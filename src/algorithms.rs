use std::{
    cmp::min,
    collections::{BTreeSet, HashSet}
};

use num::{Bounded, Zero};

use crate::types::{BlossomData, Edge, Graph, Neighborhood, RootedSubGraphForest, SquareMatrix};

pub fn floyd_warshall<W: Zero + Bounded + PartialOrd + Copy + Sized>(
    weighted_edge_list: &Vec<Edge<W>>,
    directed: bool,
) -> SquareMatrix<W> {
    let mut result =
        SquareMatrix::from_weighted_edges(weighted_edge_list, directed, W::max_value(), None);
    for i in 0..result.len {
        result.set_ix(i, i, W::zero());
    }

    //let mut pathtracker = Vec::with_capacity(result.len*result.len);
    //for i in 0..result.len {
    //    pathtracker.extend(vec![i; result.len])
    //}

    for j in 0..result.len {
        for i in 0..result.len {
            if i != j {
                for k in 0..result.len {
                    if k != j {
                        let sum = result.get_ix(i, j) + result.get_ix(j, k);
                        if result.get_ix(i, k) > sum {
                            result.set_ix(i, k, sum);
                            //pathtracker[result.len*i + k] = pathtracker[result.len*j + k];
                        }
                    }
                }
            }
        }
    }

    result
}

fn contract_blossom<W, N: Neighborhood + FromIterator<usize> + Clone>(
    graph: &Graph<W, N>,
    matching: &Vec<(usize, usize)>,
    blossom: &Vec<usize>,
) -> (Graph<W, N>, Vec<(usize, usize)>, usize, BTreeSet<usize>) {
    // assume the first vertex in the blossom is the root
    // the blossom passed to this function should always satisfy this assumption
    let blossom_root = blossom[0];

    let mut new_vertices = graph.vertices.clone();
    for v in blossom.iter() {
        if *v != blossom_root {
            new_vertices.remove(v);
        }
    }

    // construct the adjacency list of the contracted graph
    // by connecting every node that is connected to the blossom
    // directly to the root of the blossom
    let mut new_adjacencies = vec![N::new(); graph.vertices.len()];
    for i in 0..new_adjacencies.len() {
        if i == blossom_root {
            new_adjacencies[i] = graph.adjacency_list[i].clone();
            for b in blossom.iter() {
                for v in graph.adjacency_list[*b].into_iter_no_move() {
                    if !blossom.iter().any(|b| *b == v) {
                        new_adjacencies[i].insert(v);
                    }
                }
            }
        } else if !blossom.iter().any(|b| *b == i) {
            new_adjacencies[i] = graph.adjacency_list[i].clone();
            for b in blossom.iter() {
                if *b != blossom_root {
                    if new_adjacencies[i].remove(*b) {
                        new_adjacencies[i].insert(blossom_root);
                    }
                }
            }
        }
    }

    let mut new_matching = Vec::new();
    for (v, w) in matching.iter() {
        let v_in_blossom = blossom.iter().any(|u| *u == *v);
        let w_in_blossom = blossom.iter().any(|u| *u == *w);
        //////dbg!((v, w, v_in_blossom, w_in_blossom));
        if !v_in_blossom && w_in_blossom {
            new_matching.push((*v, blossom_root));
        } else if v_in_blossom && !w_in_blossom {
            new_matching.push((blossom_root, *w));
        } else if !v_in_blossom && !w_in_blossom {
            new_matching.push((*v, *w));
        }
    }

    let mut new_exposed_vertices = new_vertices.clone();
    for (v, w) in new_matching.iter() {
        new_exposed_vertices.remove(v);
        new_exposed_vertices.remove(w);
    }

    (
        Graph {
            vertices: new_vertices,
            adjacency_list: new_adjacencies,
            edges: None,
            weight_matrix: None,
            directed: graph.directed,
        },
        new_matching,
        blossom_root,
        new_exposed_vertices,
    )
}

fn expand_blossom<W, N: Neighborhood + FromIterator<usize> + Clone>(
    original_graph: &Graph<W, N>,
    original_matching: &Vec<(usize, usize)>,
    contracted_matching: &Vec<(usize, usize)>,
    path: Vec<usize>,
    blossom: &Vec<usize>,
    blossom_root: usize,
) -> Vec<usize> {
    // this function expands the path from the graph with the contracted blossom
    // to the original graph. if the path does not pass through the blossom, we
    // do not need to do anything
    if let Some(i) = path.iter().position(|v| *v == blossom_root) {
        // if the path passes through the blossom, there are two ways to expand the path
        // and we have to carefully choose the one which ensures that the resulting path
        // is still alternating
        // begin by cutting the path into the segments before and after the blossom
        // the `first_segment` is mutable because we will append the rest of the
        // edges of the augmenting path onto it
        let mut first_segment = if i > 0 {
            Vec::from(&path[..i])
        } else {
            Vec::new()
        };
        let second_segment = if i < path.len() - 1 {
            Vec::from(&path[(i + 1)..])
        } else {
            Vec::new()
        };

        // the blossom expansion can then be subdivided into two cases
        // 1) the incoming edge to the blossom is matched or the blossom
        // is the first node in the path
        // 2) the outgoing edge to the blossom is matched or the blossom
        // is the last node in the path
        let incoming_edge_matched = if i == 0 {
            true
        } else if i == path.len() - 1 {
            false
        } else {
            let incoming_vert = *first_segment.last().unwrap();
            contracted_matching
                .iter()
                .any(|e| *e == (incoming_vert, blossom_root) || *e == (blossom_root, incoming_vert))
        };

        if incoming_edge_matched {
            let outgoing_vert = second_segment[0];
            let outgoing_blossom_vert_idx =
                    match blossom.into_iter().position(|v| original_graph.adjacency_list[outgoing_vert].contains(*v)) {
                        Some(i) => i,
                        None => panic!("Original graph adjacency list did not find a blossom vertex for incoming vertex.")
                    };
            let outgoing_blossom_vert = blossom[outgoing_blossom_vert_idx];

            // absolutely obnoxious case analysis here
            for edge in original_matching.into_iter() {
                if edge.0 == outgoing_blossom_vert {
                    if edge.1 == blossom[(outgoing_blossom_vert_idx + 1) % blossom.len()] {
                        let mut intermediate_blossom_verts =
                            Vec::from(&blossom[outgoing_blossom_vert_idx..]);
                        intermediate_blossom_verts.push(blossom_root);
                        // reverse the intermediate vertices because they are starting from the "outgoing" vertex
                        first_segment.extend(intermediate_blossom_verts.iter().rev());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    } else if edge.1
                        == blossom[(outgoing_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
                    {
                        let intermediate_blossom_verts =
                            Vec::from(&blossom[0..outgoing_blossom_vert_idx + 1]);
                        first_segment.extend(intermediate_blossom_verts.iter().rev());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    }
                } else if edge.1 == outgoing_blossom_vert {
                    if edge.0 == blossom[(outgoing_blossom_vert_idx + 1) % blossom.len()] {
                        let mut intermediate_blossom_verts =
                            Vec::from(&blossom[outgoing_blossom_vert_idx..]);
                        intermediate_blossom_verts.push(blossom_root);
                        // reverse the intermediate vertices because they are starting from the "outgoing" vertex
                        first_segment.extend(intermediate_blossom_verts.iter().rev());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    } else if edge.0
                        == blossom[(outgoing_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
                    {
                        let intermediate_blossom_verts =
                            Vec::from(&blossom[0..outgoing_blossom_vert_idx + 1]);
                        first_segment.extend(intermediate_blossom_verts.iter().rev());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    }
                }
            }
            panic!("Did not find matched edge in blossom connected to outgoing vertex.")
        } else {
            let incoming_vert = *first_segment.last().unwrap();
            let incoming_blossom_vert_idx =
                    match blossom.into_iter().position(|v| original_graph.adjacency_list[incoming_vert].contains(*v)) {
                        Some(v) => v,
                        None => panic!("Original graph adjacency list did not find a blossom vertex for incoming vertex.")
                    };
            let incoming_blossom_vert = blossom[incoming_blossom_vert_idx];

            for edge in original_matching.into_iter() {
                if edge.0 == incoming_blossom_vert {
                    if edge.1 == blossom[(incoming_blossom_vert_idx + 1) % blossom.len()] {
                        let mut intermediate_blossom_verts =
                            Vec::from(&blossom[incoming_blossom_vert_idx..]);
                        intermediate_blossom_verts.push(blossom_root);
                        // reverse the intermediate vertices because they are starting from the "incoming" vertex
                        first_segment.extend(intermediate_blossom_verts.iter());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    } else if edge.1
                        == blossom[(incoming_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
                    {
                        let intermediate_blossom_verts =
                            Vec::from(&blossom[0..incoming_blossom_vert_idx + 1]);
                        first_segment.extend(intermediate_blossom_verts.iter());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    }
                } else if edge.1 == incoming_blossom_vert {
                    if edge.0 == blossom[(incoming_blossom_vert_idx + 1) % blossom.len()] {
                        let mut intermediate_blossom_verts =
                            Vec::from(&blossom[incoming_blossom_vert_idx..]);
                        intermediate_blossom_verts.push(blossom_root);
                        // reverse the intermediate vertices because they are starting from the "incoming" vertex
                        first_segment.extend(intermediate_blossom_verts.iter());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    } else if edge.0
                        == blossom[(incoming_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
                    {
                        let intermediate_blossom_verts =
                            Vec::from(&blossom[0..incoming_blossom_vert_idx + 1]);
                        first_segment.extend(intermediate_blossom_verts.iter());
                        first_segment.extend(second_segment.iter());
                        return first_segment;
                    }
                }
            }
            panic!("Did not find matched edge in blossom connected to incoming vertex.")
        }
    } else {
        path
    }
}

fn combine_root_paths(path1: &Vec<usize>, path2: &Vec<usize>) -> Vec<usize> {
    // combine two possibly overlapping paths into a single cycle
    let min_len = min(path1.len(), path2.len());
    let mut blossom: Vec<usize> = path1.into_iter().map(|n| *n).rev().collect();
    if min_len < 2 {
        return blossom;
    }
    let mut insertion_idx = 0;
    for i in (1..(min_len + 1)).rev() {
        if blossom[blossom.len() - i] == path2[i - 1] {
            insertion_idx = i;
            break;
        }
    }
    let l = blossom.len();
    for i in 0..(path2.len() - insertion_idx) {
        if i < insertion_idx {
            blossom[l + i - insertion_idx] = path2[i + insertion_idx - 1];
        } else {
            blossom.push(path2[i + insertion_idx - 1]);
        }
    }
    blossom
}

fn find_augmenting_path<W, N: Neighborhood + FromIterator<usize> + Clone>(
    graph: &Graph<W, N>,
    matching: &Vec<(usize, usize)>,
    exposed_vertices: &BTreeSet<usize>,
) -> Vec<usize> {
    // start with no marked vertices
    let mut marked_vertices = HashSet::new();
    // mark every edge that is already in the matching
    let mut marked_edges = HashSet::new();
    for edge in matching.iter() {
        marked_edges.insert(*edge);
        marked_edges.insert((edge.1, edge.0));
    }
    // roots of the subgraph forest are the exposed vertices
    // we are looking for augmenting paths, which will be paths between exposed vertices
    let mut forest: RootedSubGraphForest<N> =
        RootedSubGraphForest::new(exposed_vertices.into_iter_no_move().collect());

    while let Some(v) = forest.check_for_unmarked_vertex_w_even_dist(&marked_vertices) {
        while let Some(w) = graph.check_for_unmarked_edge(&marked_edges, v) {
            match forest.maybe_get_distance(w) {
                None => {
                    // if w is not in the forest then it is not exposed and has some edge in the matching
                    // find that edge in the matching and add it to the forest
                    let w_matched_edge = match matching.iter().find(|e| e.0 == w || e.1 == w) {
                        Some(e) => if e.0 == w {*e} else {(e.1, e.0)},
                        None => panic!("Found a vertex incident to an unmarked edge which was not in the forest and was also not in the matching.")
                    };
                    // it is important that the edges are added in this order
                    // the `add_edge` function always expects the first vertex to be in the forest already
                    forest.add_edge((v, w));
                    forest.add_edge(w_matched_edge);
                }
                Some(d) => {
                    if d % 2 == 0 {
                        // we can only find an augmenting path between v and w if the distance between them is even
                        let v_index = forest.get_vertex_index(v);
                        let w_index = forest.get_vertex_index(w);
                        if forest.root_map[v_index] != forest.root_map[w_index] {
                            // if the roots of v and w are distinct, the path between the roots is an augmenting path
                            let mut vw_path = forest.root_paths[v_index].clone();
                            vw_path.extend(forest.root_paths[w_index].iter().rev());
                            return vw_path;
                        } else {
                            // otherwise, it is an odd cycle
                            // i.e. a blossom
                            let blossom = combine_root_paths(
                                &forest.root_paths[v_index],
                                &forest.root_paths[w_index],
                            );
                            let (
                                contracted_graph,
                                contracted_matching,
                                blossom_root,
                                new_exposed_vertices,
                            ) = contract_blossom(graph, matching, &blossom);
                            let contracted_path = find_augmenting_path(
                                &contracted_graph,
                                &contracted_matching,
                                &new_exposed_vertices,
                            );
                            return expand_blossom(
                                graph,
                                matching,
                                &contracted_matching,
                                contracted_path,
                                &blossom,
                                blossom_root,
                            );
                        }
                    }
                }
            }
            marked_edges.insert((v, w));
            marked_edges.insert((w, v));
        }
        marked_vertices.insert(v);
    }
    Vec::new()
}

pub fn max_cardinality_matching<W, N: Neighborhood + FromIterator<usize> + Clone>(
    graph: &Graph<W, N>,
) -> Vec<(usize, usize)> {
    let mut matching = Vec::new();
    let mut exposed_vertices = graph.vertices.clone();

    while exposed_vertices.len() > 0 {
        let augmenting_path = find_augmenting_path(&graph, &matching, &exposed_vertices);
        if augmenting_path.len() == 0 {
            break;
        } else {
            // endpoints of the augmenting path will no longer be exposed
            let l = augmenting_path.len();
            let v = augmenting_path[0];
            let w = augmenting_path[l - 1];
            exposed_vertices.remove(&v);
            exposed_vertices.remove(&w);

            // now augment the matching, in place
            for i in 0..(augmenting_path.len() - 1) / 2 {
                if let Some(j) = matching.iter().position(|e| {
                    *e == (augmenting_path[2 * i + 1], augmenting_path[2 * i + 2])
                        || *e == (augmenting_path[2 * i + 2], augmenting_path[2 * i + 1])
                }) {
                    matching[j] = (augmenting_path[2 * i], augmenting_path[2 * i + 1]);
                } else {
                    panic!("Did not find an odd edge of the augmenting path in the matching.")
                }
            }
            matching.push((augmenting_path[l - 2], augmenting_path[l - 1]));
        }
    }
    matching
}

pub fn max_weight_matching(
    weighted_edges: &Vec<(usize, usize, f64)>,
    max_cardinality: bool,
) -> Vec<i64> {
    // not modified during algorithm iteration
    let mut n_vertices = 0;
    let mut max_weight = f64::min_value();
    for &edge in weighted_edges.iter() {
        if edge.0 == edge.1 {
            panic!("Invalid edge {:?}", edge);
        }
        if edge.0 >= n_vertices {
            n_vertices = edge.0 + 1;
        }
        if edge.1 >= n_vertices {
            n_vertices = edge.1 + 1;
        }
        if edge.2 > max_weight {
            max_weight = edge.2;
        }
    }
    let mut neighborhood_endpoints = vec![Vec::new(); n_vertices];
    for (i, edge) in weighted_edges.iter().enumerate() {
        neighborhood_endpoints[edge.0].push(2 * i + 1);
        neighborhood_endpoints[edge.1].push(2 * i);
    }
    let mut endpoints = Vec::with_capacity(2 * weighted_edges.len());
    for edge in weighted_edges.iter() {
        endpoints.push(edge.0);
        endpoints.push(edge.1);
    }

    // modified during algorithm iteration
    let mut blossom_data = BlossomData::new(n_vertices, weighted_edges.len(), max_weight);
    let mut matching = vec![-1i64; n_vertices];

    for _stage in 0..n_vertices {
        //println!("STAGE: {} MATCHING: {:?}", _stage, matching);

        blossom_data.clear();
        let mut stack = Vec::with_capacity(n_vertices);
        //////dbg!(blossom_data.allowed_edge.clone());

        // prepare to root a search tree at every unmatched vertex
        // which is not contained in a contracted blossom
        for v in 0..n_vertices {
            let id = blossom_data.blossom_id[v] as usize;
            let label = blossom_data.blossom_labels[id];
            if matching[v] == -1 && label == 0 {
                blossom_data.assign_label(&mut stack, v, 1, -1, &endpoints, &matching);
            }
        }
        //////dbg!(blossom_data.allowed_edge.clone());

        let mut augmented = false;
        loop {
            //println!("SUBSTAGE");
            while stack.len() > 0 && !augmented {
                //println!("STACK {:?}", stack);
                let v = stack.pop().unwrap();
                //println!("POP: {}", v);
                println!("LABELS {:?}", blossom_data.blossom_labels);

                assert_eq!(blossom_data.blossom_labels[blossom_data.blossom_id[v]], 1);

                ////dbg!(neighborhood_endpoints[v].clone());
                for &p in neighborhood_endpoints[v].iter() {
                    let edge_idx = p / 2;
                    let w = endpoints[p];
                    ////println!("p v w blossom_id[v] blossom_id[w] {} {} {} {} {}", p, v, w, blossom_data.blossom_id[v], blossom_data.blossom_id[w]);
                    if blossom_data.blossom_id[v] == blossom_data.blossom_id[w] {
                        continue;
                    }

                    let mut edge_slack = 0.0;
                    //dbg!(p, edge_idx);
                    if !blossom_data.allowed_edge[edge_idx] {
                        edge_slack = blossom_data.slack(edge_idx, weighted_edges);
                        //println!("edge_slack {}", edge_slack);
                        if edge_slack <= 1e-12 {
                            blossom_data.allowed_edge[edge_idx] = true;
                        }
                    }

                    let label = blossom_data.blossom_labels[blossom_data.blossom_id[w]];
                    //println!("v w allowed_edge {} {} {}", v, w, blossom_data.allowed_edge[edge_idx]);
                    //////dbg!(label);
                    if blossom_data.allowed_edge[edge_idx] {
                        if label == 0 {
                            blossom_data.assign_label(&mut stack, w, 2, (p ^ 1) as i64, &endpoints, &matching);
                        } else if label == 1 {
                            let root = blossom_data
                                .check_for_blossom_or_augmenting_path(v, w, &endpoints, &matching);
                            match root {
                                None => {
                                    //println!("MATCHING {:?}", matching);
                                    //println!("{}", edge_idx);
                                    blossom_data.augment_matching(
                                        edge_idx,
                                        weighted_edges,
                                        &endpoints,
                                        &mut matching,
                                    );
                                    augmented = true;
                                    //println!("MATCHING {:?}", matching);
                                    break;
                                }
                                Some(r) => {
                                    blossom_data.add_blossom(
                                        &mut stack,
                                        r,
                                        edge_idx,
                                        &matching,
                                        &endpoints,
                                        &neighborhood_endpoints,
                                        &weighted_edges,
                                    );
                                }
                            }
                        } else if blossom_data.blossom_labels[w] == 0 {
                            assert_eq!(label, 2);
                            blossom_data.blossom_labels[w] = 2;
                            blossom_data.label_endpoints[w] = (p ^ 1) as i64;
                        }
                    } else if label == 1 {
                        let b = blossom_data.blossom_id[v];
                        //////dbg!(blossom_data.best_edge[b]);
                        //////dbg!((edge_slack, blossom_data.slack(blossom_data.best_edge[b] as usize, &edges, &weight_matrix)));
                        if blossom_data.best_edge[b] < 0
                            || edge_slack
                                < blossom_data.slack(
                                    blossom_data.best_edge[b] as usize,
                                    weighted_edges
                                )
                        {
                            blossom_data.best_edge[b] = edge_idx as i64;
                            //////dbg!(b, blossom_data.best_edge[b]);
                        }
                    } else if label == 0 {
                        if blossom_data.best_edge[w] < 0
                            || edge_slack
                                < blossom_data.slack(
                                    blossom_data.best_edge[w] as usize,
                                    weighted_edges
                                )
                        {
                            blossom_data.best_edge[w] = edge_idx as i64;
                        }
                    }
                }
            }
            if augmented {
                //println!("AUGMENTED");
                break;
            }

            let (delta, best_edge, update_blossom) = blossom_data.determine_delta_and_update_mode(
                weighted_edges,
                max_cardinality,
            );
            blossom_data.update_dual_soln(delta);
            let must_break = blossom_data.update_blossom_structure(
                &mut stack,
                best_edge,
                update_blossom,
                weighted_edges,
                &endpoints,
                &matching,
            );
            if must_break {
                break;
            }
        }

        // no more augmenting path can be found
        if !augmented {
            break;
        }

        // end of this stage, expand all S-blossoms which have dual_soln = 0
        for b in n_vertices..2 * n_vertices {
            if blossom_data.s_blossom_is_tight(b) {
                blossom_data.expand_blossom(&mut stack, b, true, &endpoints, &matching);
            }
        }
    }

    // as we are done traversing the graph, we transform
    // endpoint indices back into actual vertex indicies
    for v in 0..n_vertices {
        if matching[v] >= 0 {
            matching[v] = endpoints[matching[v] as usize] as i64;
        }
    }

    matching
}
