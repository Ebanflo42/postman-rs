use std::{
    cmp::min,
    collections::{BTreeSet, HashSet},
};

use num::{Bounded, Zero};

use crate::types::{Edge, Graph, Neighborhood, RootedSubGraphForest, SquareMatrix};

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
        dbg!((v, w, v_in_blossom, w_in_blossom));
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
    for i in (1..(min_len+1)).rev() {
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

pub fn edmonds_max_cardinality_matching<W, N: Neighborhood + FromIterator<usize> + Clone>(
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
