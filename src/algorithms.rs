use std::path;

use num::{Bounded, Zero, One};
use ndarray::Array2;

use crate::types::{BlossomData};

//*
pub fn floyd_warshall<W: Zero + One + Bounded + PartialOrd + Copy + Sized>(
    weighted_edges: &Vec<(usize, usize, W)>,
    directed: bool,
) -> (Array2<W>, Array2<i64>) {
    let mut n_vertices = 0;
    for &edge in weighted_edges.iter() {
        if n_vertices < edge.0 {
            n_vertices = edge.0;
        }
        if n_vertices < edge.1 {
            n_vertices = edge.1;
        }
    }

    let mut distances = Array2::from_shape_fn((n_vertices, n_vertices), |_| W::max_value());
    let mut pathtracker = Array2::from_shape_fn((n_vertices, n_vertices), |_| -1);
    for &edge in weighted_edges.iter() {
        distances[[edge.0, edge.1]] = edge.2;
        pathtracker[[edge.0, edge.1]] = edge.0 as i64;
        if !directed {
            distances[[edge.1, edge.0]] = edge.2;
            pathtracker[[edge.1, edge.0]] = edge.0 as i64;
        }
    }
    for j in 0..n_vertices {
        distances[[j, j]] = W::zero();
        pathtracker[[j, j]] = j as i64;
    }

    for j in 0..n_vertices {
        for i in 0..n_vertices {
            if i != j {
                for k in 0..n_vertices {
                    if k != j {
                        let sum = distances[[i, j]] + distances[[j, k]];
                        if distances[[i, k]] > sum {
                            distances[[i, k]] = sum;
                            pathtracker[[i, k]] = pathtracker[[j, k]];
                        }
                    }
                }
            }
        }
    }

    (distances, pathtracker)
}
//*/

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
    max_weight = f64::max(0.0, max_weight);
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

        blossom_data.clear();
        let mut stack = Vec::with_capacity(n_vertices);

        // prepare to root a search tree at every unmatched vertex
        // which is not contained in a contracted blossom
        for v in 0..n_vertices {
            let id = blossom_data.blossom_id[v] as usize;
            let label = blossom_data.blossom_labels[id];
            if matching[v] == -1 && label == 0 {
                blossom_data.assign_label(&mut stack, v, 1, -1, &endpoints, &matching);
            }
        }

        let mut augmented = false;
        let mut count = 0;
        loop {
            count += 1;
            if count > 10 {
                break;
            }
            while stack.len() > 0 && !augmented {
                let v = stack.pop().unwrap();

                //assert_eq!(blossom_data.blossom_labels[blossom_data.blossom_id[v]], 1);
                for &p in neighborhood_endpoints[v].iter() {
                    let edge_idx = p / 2;
                    let w = endpoints[p];
                    if blossom_data.blossom_id[v] == blossom_data.blossom_id[w] {
                        continue;
                    }

                    let mut edge_slack = 0.0;
                    if !blossom_data.allowed_edge[edge_idx] {
                        edge_slack = blossom_data.slack(edge_idx, weighted_edges);
                        if edge_slack <= 1e-12 {
                            blossom_data.allowed_edge[edge_idx] = true;
                        }
                    }

                    let label = blossom_data.blossom_labels[blossom_data.blossom_id[w]];
                    if blossom_data.allowed_edge[edge_idx] {
                        if label == 0 {
                            blossom_data.assign_label(&mut stack, w, 2, (p ^ 1) as i64, &endpoints, &matching);
                        } else if label == 1 {
                            let root = blossom_data
                                .check_for_blossom_or_augmenting_path(v, w, &endpoints, &matching);
                            match root {
                                None => {
                                    blossom_data.augment_matching(
                                        edge_idx,
                                        weighted_edges,
                                        &endpoints,
                                        &mut matching,
                                    );
                                    augmented = true;
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
                        if blossom_data.best_edge[b] < 0
                            || edge_slack
                                < blossom_data.slack(
                                    blossom_data.best_edge[b] as usize,
                                    weighted_edges
                                )
                        {
                            blossom_data.best_edge[b] = edge_idx as i64;
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

pub fn min_weight_max_cardinality_matching(
    weighted_edges: &Vec<(usize, usize, f64)>,
) -> Vec<i64> {
    let weighted_edges = weighted_edges.iter().map(|&(i, j, w)| (i, j, -w)).collect();
    max_weight_matching(&weighted_edges, true)
}
