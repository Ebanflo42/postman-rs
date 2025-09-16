use num::{Bounded, Zero};

use crate::types::{BlossomData};

/*
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
*/
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
            dbg!("LOOP");
            count += 1;
            if count > 10 {
                break;
            }
            while stack.len() > 0 && !augmented {
                let v = stack.pop().unwrap();

                assert_eq!(blossom_data.blossom_labels[blossom_data.blossom_id[v]], 1);
                println!("{:?}", stack);
                for &p in neighborhood_endpoints[v].iter() {
                    let edge_idx = p / 2;
                    let w = endpoints[p];
                    if blossom_data.blossom_id[v] == blossom_data.blossom_id[w] {
                        continue;
                    }

                    let mut edge_slack = 0.0;
                    if !blossom_data.allowed_edge[edge_idx] {
                        edge_slack = blossom_data.slack(edge_idx, weighted_edges);
                        //println!("edge_slack {}", edge_slack);
                        if edge_slack <= 1e-12 {
                            blossom_data.allowed_edge[edge_idx] = true;
                        }
                    }

                    let label = blossom_data.blossom_labels[blossom_data.blossom_id[w]];
                    //println!("v w allow_edge {} {} {}", v, w, blossom_data.allowed_edge[edge_idx]);
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
                println!("AUGMENTED");
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
                println!("deltatype 1");
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
        loop {
            while stack.len() > 0 && !augmented {
                let v = stack.pop().unwrap();

                assert_eq!(blossom_data.blossom_labels[blossom_data.blossom_id[v]], 1);

                for &p in neighborhood_endpoints[v].iter() {
                    let edge_idx = p / 2;
                    let w = endpoints[p];
                    if blossom_data.blossom_id[v] == blossom_data.blossom_id[w] {
                        continue;
                    }

                    let mut edge_slack = 0.0;
                    if !blossom_data.allowed_edge[edge_idx] {
                        edge_slack = blossom_data.slack(edge_idx, weighted_edges);
                        //println!("edge_slack {}", edge_slack);
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
                true,
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
