use core::fmt::Display;
use std::{collections::{BTreeSet, HashSet}, marker::PhantomData, ops::Deref};

use num::{Bounded, PrimInt, Zero};

use crate::types::{Edge, Graph, Neighborhood, RootedSubGraphForest, SquareMatrix};

pub fn floyd_warshall<W: Bounded + Zero + Ord + Copy + Deref<Target = W>>(
    weighted_edge_list: Vec<Edge<W>>,
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
                    let sum = result.get_ix(i, j) + result.get_ix(j, k);
                    if k != j && result.get_ix(i, k) > sum {
                        result.set_ix(i, k, sum);
                        //pathtracker[result.len*i + k] = pathtracker[result.len*j + k];
                    }
                }
            }
        }
    }

    result
}

fn contract_blossom<
    W,
    N: Neighborhood + FromIterator<usize> + IntoIterator<Item = usize> + Clone,
>(
    graph: &Graph<W, N>,
    matching: &Vec<(usize, usize)>,
    blossom: &Vec<usize>,
) -> (Graph<W, N>, Vec<(usize, usize)>, usize) {
    // assume the first vertex in the blossom is the root
    // the blossom passed to this function should always satisfy this assumption
    let blossom_root = blossom[0];

    let mut new_vertices = graph.vertices.clone();
    for v in blossom.iter() {
        if *v != blossom_root {
            new_vertices.remove(v);
        }
    }

    let mut new_adjacencies = vec![N::new(); graph.vertices.len()];
    for i in 0..new_adjacencies.len() {
        if i == blossom_root {
            new_adjacencies[i] = graph.adjacency_list[i].clone();
            for b in blossom.iter() {
                for v in graph.adjacency_list[*b].into_iter() {
                    new_adjacencies[i].insert(v);
                }
            }
        } else if !blossom.iter().any(|b| *b == i) {
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
        if !v_in_blossom && w_in_blossom {
            new_matching.push((*v, blossom_root));
        } else if v_in_blossom && !w_in_blossom {
            new_matching.push((blossom_root, *w));
        } else if !v_in_blossom && !w_in_blossom {
            new_matching.push((*v, *w));
        }
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
    )
}

fn expand_blossom<W, N: Neighborhood + FromIterator<usize> + IntoIterator<Item = usize> + Clone>(
    original_graph: &Graph<W, N>,
    original_matching: &Vec<(usize, usize)>,
    contracted_matching: &Vec<(usize, usize)>,
    path: Vec<usize>,
    blossom: &Vec<usize>,
    blossom_root: usize,
) -> Vec<usize> {
    if let Some(i) = path.iter().position(|v| *v == blossom_root) {
        // if the path passes through the blossom, there are two ways to expand the path
        // and we have to carefully choose the one which ensures that the resulting path
        // is still alternating
        // begin by cutting the path into the segments before and after the blossom
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
        let incoming_vert = *first_segment.last().unwrap();
        let outgoing_vert = second_segment[0];

        // determine if the edge going into the blossom is in the matching
        // or if the blossom is the first node in the path
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

fn find_augmenting_path<
    W,
    N: Neighborhood + FromIterator<usize> + IntoIterator<Item = usize> + Clone,
>(
    graph: &Graph<W, N>,
    matching: &Vec<(usize, usize)>,
    exposed_vertices: &BTreeSet<usize>,
) -> Vec<usize> {
    let mut marked_vertices = HashSet::new();
    let mut marked_edges = HashSet::from_iter(matching.iter().map(|x| *x));
    let mut forest: RootedSubGraphForest<N> = RootedSubGraphForest::new(exposed_vertices.iter().map(|x| *x).collect());

    while let Some(v) = forest.check_for_unmarked_vertex_w_even_dist(&marked_vertices) {
        while let Some(w) = forest.check_for_unmarked_edge(&marked_edges, v) {
            match forest.maybe_get_distance(w) {
                None => {
                    let w_matched_edge = match matching.iter().position(|e| e.0 == w || e.1 == w) {
                        Some(e) => matching[e],
                        None => panic!("Found a vertex incident to an unmarked edge which was not in the forest and was also not in the matching.")
                    };
                    let x = if w_matched_edge.0 == w {
                        w_matched_edge.1
                    } else {
                        w_matched_edge.0
                    };
                    forest.add_edge((v, w));
                    forest.add_edge((w, x));
                }
                Some(d) => {
                    if d % 2 == 0 {
                        let v_index = forest.get_vertex_index(v);
                        let w_index = forest.get_vertex_index(w);
                        let mut vw_path = forest.root_paths[v_index].clone();
                        vw_path.extend(forest.root_paths[w_index].iter().rev());
                        if forest.root_map[v_index] != forest.root_map[w_index] {
                            return vw_path;
                        } else {
                            let (contracted_graph, contracted_matching, blossom_root) =
                                contract_blossom(graph, matching, &vw_path);
                            let contracted_path = find_augmenting_path(
                                &contracted_graph,
                                &contracted_matching,
                                exposed_vertices,
                            );
                            return expand_blossom(
                                graph,
                                matching,
                                &contracted_matching,
                                contracted_path,
                                &vw_path,
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

pub fn edmonds_max_cardinality_matching<
    W,
    N: Neighborhood + FromIterator<usize> + IntoIterator<Item = usize> + Clone,
>(
    graph: &Graph<W, N>,
) -> Vec<(usize, usize)> {
    let mut matching = Vec::new();
    let mut exposed_vertices = graph.vertices.clone();
    loop {
        let augmenting_path = find_augmenting_path(&graph, &matching, &exposed_vertices);
        if augmenting_path.len() == 0 {
            break;
        } else {
            // endpoints of the augmenting path will no longer be exposed
            let l = augmenting_path.len();
            let v = augmenting_path[0];
            let w = augmenting_path[l - 1];
            if let Some(i) = exposed_vertices.iter().position(|u| *u == v) {
                exposed_vertices.remove(&i);
            } else {
                panic!("First vertex of augmenting path was not an exposed vertex!");
            }
            if let Some(j) = exposed_vertices.iter().position(|u| *u == w) {
                exposed_vertices.remove(&j);
            } else {
                panic!("Last vertex of augmenting path was not an exposed vertex!");
            }

            // now augment the matching, in place
            for i in 0..(augmenting_path.len() - 1) / 2 {
                if let Some(j) = matching.iter().position(|e| {
                    *e == (augmenting_path[2 * i + 1], augmenting_path[2 * i + 2])
                        || *e == (augmenting_path[2 * i + 2], augmenting_path[2 * i + 1])
                }) {
                    matching[j] = (augmenting_path[2 * i], augmenting_path[2 * i + 1]);
                } else {
                    panic!("Did not find an edge of the augmenting path in the matching.")
                }
            }
            matching.push((augmenting_path[l - 2], augmenting_path[l - 1]));
        }
    }
    matching
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floyd_warshall1() {
        let weighted_edges = vec![
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 2.0),
            (3, 4, 1.0),
            (4, 0, 1.0),
            (0, 3, 1.0),
        ];
        let distance_matrix = floyd_warshall(weighted_edges, false);
        println!("{}", distance_matrix.to_string());
        assert_eq!(
            distance_matrix.data,
            vec![
                0.0, 1.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 2.0, 3.0, 1.0,
                2.0, 2.0, 0.0, 1.0, 1.0, 2.0, 3.0, 1.0, 0.0
            ]
        );
    }
}
