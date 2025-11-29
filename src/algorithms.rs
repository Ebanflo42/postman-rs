//#![postman = "doc"]

use std::collections::{BTreeSet, HashSet};
use std::fmt::Debug;

use ndarray::Array2;
use num::{Bounded, Num, NumCast};

use crate::blossom_data::{BlossomData, BlossomLabel};

//*
pub fn floyd_warshall<W: Bounded + PartialOrd + Copy + Sized + Num + NumCast>(
    weighted_edges: &Vec<(usize, usize, W)>,
    directed: bool,
) -> (Array2<W>, Array2<i64>) {
    /// The Floyd-Warshall algorithm.
    ///
    /// Takes as input a list of edges of the form (u, v, w) where u and v are vertices
    /// and w is a weight. Assumes vertices from 0 to the maximum vertex supplied in the
    /// edges exist.
    ///
    /// Input `directed` specifies whether or not the edges are considered as directed or
    /// undirected.
    ///
    /// The user is responsible for ensuring duplicate edges are not given.
    ///
    /// Returns `(distances, pathtracker)`, where `distances[[u, v]]` is the length of the
    /// shortest path between vertices u and v and `pathtracker[[u, v]]` is the first step
    /// in the shortest path from from v to u.
    let mut n_vertices = 0;
    for &edge in weighted_edges.iter() {
        if n_vertices <= edge.0 {
            n_vertices = edge.0 + 1;
        }
        if n_vertices <= edge.1 {
            n_vertices = edge.1 + 1;
        }
    }

    let mut distances = Array2::from_shape_fn((n_vertices, n_vertices), |_| {
        W::max_value() / W::from(2).unwrap()
    });
    let mut pathtracker = Array2::from_shape_fn((n_vertices, n_vertices), |_| -1);
    for &edge in weighted_edges.iter() {
        distances[[edge.0, edge.1]] = edge.2;
        pathtracker[[edge.0, edge.1]] = edge.0 as i64;
        if !directed {
            distances[[edge.1, edge.0]] = edge.2;
            pathtracker[[edge.1, edge.0]] = edge.1 as i64;
        }
    }
    for j in 0..n_vertices {
        distances[[j, j]] = W::zero();
        pathtracker[[j, j]] = j as i64;
    }

    for k in 0..n_vertices {
        for i in 0..n_vertices {
            if i != k {
                for j in 0..n_vertices {
                    if j != i && j != k {
                        let sum = distances[[i, k]] + distances[[k, j]];
                        if distances[[i, j]] > sum {
                            distances[[i, j]] = sum;
                            pathtracker[[i, j]] = pathtracker[[k, j]];
                        }
                    }
                }
            }
        }
    }

    (distances, pathtracker)
}

fn floyd_warshall2<W: Bounded + PartialOrd + Copy + Sized + Num + NumCast>(
    weighted_edges: &Vec<(usize, usize, W)>,
    n_vertices: usize,
    directed: bool,
) -> (Array2<W>, Array2<i64>) {
    /// Same as above but the number of vertices has already been computed.
    let mut distances = Array2::from_shape_fn((n_vertices, n_vertices), |_| {
        W::max_value() / W::from(2).unwrap()
    });
    let mut pathtracker = Array2::from_shape_fn((n_vertices, n_vertices), |_| -1);
    for &edge in weighted_edges.iter() {
        distances[[edge.0, edge.1]] = edge.2;
        pathtracker[[edge.0, edge.1]] = edge.0 as i64;
        if !directed {
            distances[[edge.1, edge.0]] = edge.2;
            pathtracker[[edge.1, edge.0]] = edge.1 as i64;
        }
    }
    for j in 0..n_vertices {
        distances[[j, j]] = W::zero();
        pathtracker[[j, j]] = j as i64;
    }

    for k in 0..n_vertices {
        for i in 0..n_vertices {
            if i != k {
                for j in 0..n_vertices {
                    if j != i && j != k {
                        let sum = distances[[i, k]] + distances[[k, j]];
                        if distances[[i, j]] > sum {
                            distances[[i, j]] = sum;
                            pathtracker[[i, j]] = pathtracker[[k, j]];
                        }
                    }
                }
            }
        }
    }

    (distances, pathtracker)
}

pub fn max_weight_matching<W: Bounded + PartialOrd + Copy + Sized + Num + NumCast + Into<f64>>(
    weighted_edges: &Vec<(usize, usize, W)>,
    max_cardinality: bool,
) -> Vec<i64> {
    /// Edmond's blossom algorithm for maximum weight matchings.
    ///
    /// Takes as input a list of edges of the form (u, v, w) where u and v are vertices
    /// and w is a weight. Assumes vertices from 0 to the maximum vertex supplied in the
    /// edges exist.
    ///
    /// Input `max_cardinality` should be `true` if only maximum cardinality matchings should
    /// be considered as solutions.
    ///
    /// Returns a list of "mates" to each given vertex. Vertex u is matched to mate[u] if
    /// mate[u] >= 0 and u is unmatched if mate[u] = -1.

    // not modified during algorithm iteration
    let mut n_vertices = 0;
    let mut max_weight = W::min_value();
    for &edge in weighted_edges.iter() {
        if edge.0 == edge.1 {
            panic!("Invalid edge {:?}", (edge.0, edge.1));
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
    // the max_weight is used to initialize dual variables, so clip it to zero
    // if it is negative, to ensure feasability
    max_weight = if max_weight < W::zero() {
        W::zero()
    } else {
        max_weight
    };
    // this table stores indices of edges incident to a given vertex
    let mut neighborhood_endpoints = vec![Vec::new(); n_vertices];
    for (i, edge) in weighted_edges.iter().enumerate() {
        neighborhood_endpoints[edge.0].push(2 * i + 1);
        neighborhood_endpoints[edge.1].push(2 * i);
    }
    // this table stores the endpoints of given edges
    let mut endpoints = Vec::with_capacity(2 * weighted_edges.len());
    for edge in weighted_edges.iter() {
        endpoints.push(edge.0);
        endpoints.push(edge.1);
    }

    // modified during algorithm iteration
    // this data structure holds the entire state of the "blossom algorithm"
    // this essentially amounts to a forest whose leaves are vertices and nodes share
    // a parent if they are found to coexist in an odd cycle (blossom)
    let mut blossom_data = BlossomData::new(n_vertices, weighted_edges.len(), max_weight);
    // this array represents the matching, -1 for unmatched vertices,
    // non-negative integers indicate the index of the matched vertex
    let mut matching = vec![-1i64; n_vertices];

    // at each stage, either one edge is added to the matching
    // or it is determined that no max cardinality matching exists
    // or it is determined that we found the maximum weight matching
    // therefore, the algorithm takes at most n_vertices/2 stages
    for _stage in 0..n_vertices {
        // re-initialize the blossom data without reallocating
        blossom_data.clear();
        // this stack will contain the different roots of our search forest
        let mut stack = Vec::with_capacity(n_vertices);

        for v in 0..n_vertices {
            // get the top-level blossom containing vertex v
            let id = blossom_data.blossom_id[v] as usize;
            // check the label of this blossom
            let label = blossom_data.blossom_labels[id];
            // if the vertex is unmatched and the blossom is unlabeled (unmatched)
            // assign the label "inner" to the vertex and append the children of the
            // blossom to the stack
            if matching[v] == -1 && label == BlossomLabel::Unlabeled {
                blossom_data.assign_label(&mut stack, v, BlossomLabel::Inner, -1, &endpoints, &matching);
            }
        }

        // if we augment the matching during this stage, we will proceed to the next stage
        let mut augmented = false;
        loop {
            while stack.len() > 0 && !augmented {
                let v = stack.pop().unwrap();

                //assert_eq!(blossom_data.blossom_labels[blossom_data.blossom_id[v]], 1);
                // check all the edges incident to the current vertex v
                for &p in neighborhood_endpoints[v].iter() {
                    let edge_idx = p / 2;
                    let w = endpoints[p];
                    // if the two vertices on this edge belong to the same top-level blossom,
                    // skip it
                    if blossom_data.blossom_id[v] == blossom_data.blossom_id[w] {
                        continue;
                    }

                    // check if the edge is "allowed"
                    // meaning whether or not its "slack" is zero
                    let mut edge_slack = W::zero();
                    if !blossom_data.allowed_edge[edge_idx] {
                        edge_slack = blossom_data.slack(edge_idx, weighted_edges);
                        if edge_slack.into() <= 1e-12f64 {
                            blossom_data.allowed_edge[edge_idx] = true;
                        }
                    }

                    // check the label of the top-level blossom containing vertex w
                    let label = blossom_data.blossom_labels[blossom_data.blossom_id[w]];
                    if blossom_data.allowed_edge[edge_idx] {
                        // if w hasn't been labeled, assign it label "outer"
                        // recall that the roots of our search forest are all "inner"
                        if label == BlossomLabel::Unlabeled {
                            blossom_data.assign_label(
                                &mut stack,
                                w,
                                BlossomLabel::Outer,
                                (p ^ 1) as i64,
                                &endpoints,
                                &matching,
                            );
                        } else if label == BlossomLabel::Inner {
                            // if vertex w is already labeled as "inner"
                            // we should look for either an augmenting path,
                            // which is an alternating path which will get us a better matching,
                            // or a blossom, which is an alternating path that is an odd cycle
                            let root = blossom_data
                                .check_for_blossom_or_augmenting_path(v, w, &endpoints, &matching);
                            // root is None if we found an augmented path
                            match root {
                                None => {
                                    // root is None if we found an augmenting path
                                    // then we are done with this stage
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
                                    // root is a vertex index if we found a blossom
                                    // the root of the blossom is the vertex which is
                                    // already matched to some external vertex
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
                        } else if blossom_data.blossom_labels[w] == BlossomLabel::Unlabeled {
                            // if we find that vertex w is unlabeled, we should have already
                            // matched the top-level blossom containing it
                            //assert_eq!(label, BlossomLabel::Outer);
                            // in that case, label w as outer
                            blossom_data.blossom_labels[w] = BlossomLabel::Outer;
                            blossom_data.label_endpoints[w] = (p ^ 1) as i64;
                        }
                    } else if label == BlossomLabel::Inner {
                        // if the blossom containing vertex w is already labeled as inner
                        // we might need to update the "best edge" coming from top-level
                        // blossom b containing v
                        // the "best edges" are simply the least-slack edges
                        // coming out of a vertex/blossom and it is important to
                        // record them in order to quickly compute dual updates
                        let b = blossom_data.blossom_id[v];
                        if blossom_data.best_edge[b] < 0
                            || edge_slack
                                < blossom_data
                                    .slack(blossom_data.best_edge[b] as usize, weighted_edges)
                        {
                            blossom_data.best_edge[b] = edge_idx as i64;
                        }
                    } else if label == BlossomLabel::Unlabeled {
                        // additionally, if the blossom containing w is unlabeled
                        // we may need to update the "best edge" coming from w
                        if blossom_data.best_edge[w] < 0
                            || edge_slack
                                < blossom_data
                                    .slack(blossom_data.best_edge[w] as usize, weighted_edges)
                        {
                            blossom_data.best_edge[w] = edge_idx as i64;
                        }
                    }
                }
            }
            // if we already augmented the matching, we are done with this stage
            if augmented {
                break;
            }

            let (delta, new_tight_edge, new_tight_blossom) =
                blossom_data.determine_delta_and_update_mode(weighted_edges, max_cardinality);
            blossom_data.update_dual_soln(delta);
            let optimal_solution_found = blossom_data.update_blossom_structure(
                &mut stack,
                new_tight_edge,
                new_tight_blossom,
                weighted_edges,
                &endpoints,
                &matching,
            );
            if optimal_solution_found {
                break;
            }
        }

        // if we end the stage without augmenting the matching then
        // no more augmenting path can be found and the algorithm is done
        if !augmented {
            break;
        }

        // end of this stage, expand all S-blossoms which have dual_soln = 0
        for b in n_vertices..2 * n_vertices {
            if blossom_data.inner_blossom_is_tight(b) {
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

pub fn min_weight_max_cardinality_matching<
    W: Bounded + PartialOrd + Copy + Sized + Num + NumCast + Into<f64>,
>(
    weighted_edges: &Vec<(usize, usize, W)>,
) -> Vec<i64> {
    /// Call Edmond's blossom algorithm above, but to get a minimum weight
    /// maximum cardinality matching.
    ///
    /// Returns a list of "mates" to each given vertex. Vertex u is matched to mate[u] if
    /// mate[u] >= 0 and u is unmatched if mate[u] = -1.
    let weighted_edges = weighted_edges
        .iter()
        .map(|&(i, j, w)| (i, j, W::from(-1).unwrap() * w))
        .collect();
    max_weight_matching(&weighted_edges, true)
}

pub fn min_weight_t_join<W: Bounded + PartialOrd + Copy + Sized + Num + NumCast + Into<f64>>(
    weighted_edges: &Vec<(usize, usize, W)>,
    t: &Vec<usize>,
) -> BTreeSet<(usize, usize)> {
    /// Solves the minimum weight T-join problem for the weighted graph specified by
    /// `weighted_edges` and the vertex set specified by `t`.
    ///
    /// Throws an error if `t` has odd length, and assumes there are no duplicate
    /// vertices in `t`.
    ///
    /// This function returns a `BTreeSet` of edges of the form (u, v) where u and v
    /// are `usize`s pointing to specific vertices. BTreeSets are useful here for the
    /// balance between efficient iteration over the elements and efficient checking
    /// for membership of an element.

    if t.len() % 2 != 0 {
        panic!("Vertex set `t` must have even cardinality!");
    }

    let mut n_vertices = 0;
    for &edge in weighted_edges.iter() {
        if n_vertices <= edge.0 {
            n_vertices = edge.0 + 1;
        }
        if n_vertices <= edge.1 {
            n_vertices = edge.1 + 1;
        }
    }

    let (distances, pathtracker) = floyd_warshall2(weighted_edges, n_vertices, false);
    let construct_path = |i: usize, j: usize| {
        let mut result = vec![(pathtracker[[i, j]] as usize, j)];
        while result.last().unwrap().0 != i {
            let next_step = pathtracker[[i, result.last().unwrap().0]];
            if next_step < 0 {
                panic!("{} and {} are disconnected!", i, j);
            }
            result.push((next_step as usize, result.last().unwrap().0));
        }
        result
    };

    let mut metric_closure_edges = Vec::new();
    for i in 0..t.len() {
        for j in 0..i {
            metric_closure_edges.push((i, j, distances[[t[i], t[j]]]));
        }
    }

    let metric_closure_matching = min_weight_max_cardinality_matching(&metric_closure_edges);

    let mut covered_indices = HashSet::new();
    let mut result_edges = BTreeSet::new();
    for (i, &j) in metric_closure_matching.iter().enumerate() {
        if j < 0 {
            panic!("Did not find a perfect matching on the metric closure of `t`");
        }
        if covered_indices.contains(&i) || covered_indices.contains(&(j as usize)) {
            continue;
        }

        let path = construct_path(t[i], t[j as usize]);
        // symmetric difference accounting for the fact that
        // edges are unordered tuples
        for &edge in path.iter() {
            if result_edges.contains(&edge) {
                result_edges.remove(&edge);
            } else if result_edges.contains(&(edge.1, edge.0)) {
                result_edges.remove(&(edge.1, edge.0));
            } else {
                result_edges.insert(edge);
            }
        }

        covered_indices.insert(i);
        covered_indices.insert(j as usize);
    }

    result_edges
}

fn min_weight_t_join2<
    W: Debug + Bounded + PartialOrd + Copy + Sized + Num + NumCast + Into<f64>,
>(
    weighted_edges: &Vec<(usize, usize, W)>,
    t: &Vec<usize>,
    n_vertices: usize,
) -> BTreeSet<(usize, usize)> {
    /// Same as above, but the number of vertices is supplied as an argument.

    if t.len() % 2 != 0 {
        panic!("Vertex set `t` must have even cardinality!");
    }

    let (distances, pathtracker) = floyd_warshall2(weighted_edges, n_vertices, false);
    let construct_path = |i: usize, j: usize| {
        let mut result = vec![(pathtracker[[i, j]] as usize, j)];
        while result.last().unwrap().0 != i {
            let next_step = pathtracker[[i, result.last().unwrap().0]];
            if next_step < 0 {
                panic!("{} and {} are disconnected!", i, j);
            }
            result.push((next_step as usize, result.last().unwrap().0));
        }
        result
    };

    let mut metric_closure_edges = Vec::new();
    for i in 0..t.len() {
        for j in 0..i {
            metric_closure_edges.push((i, j, distances[[t[i], t[j]]]));
        }
    }

    let metric_closure_matching = min_weight_max_cardinality_matching(&metric_closure_edges);

    let mut covered_indices = HashSet::new();
    let mut result_edges = BTreeSet::new();
    for (i, &j) in metric_closure_matching.iter().enumerate() {
        if j < 0 {
            panic!("Did not find a perfect matching on the metric closure of `t`");
        }
        if covered_indices.contains(&i) || covered_indices.contains(&(j as usize)) {
            continue;
        }

        let path = construct_path(t[i], t[j as usize]);
        // symmetric difference accounting for the fact that
        // edges are unordered tuples
        for &edge in path.iter() {
            if result_edges.contains(&edge) {
                result_edges.remove(&edge);
            } else if result_edges.contains(&(edge.1, edge.0)) {
                result_edges.remove(&(edge.1, edge.0));
            } else {
                result_edges.insert(edge);
            }
        }

        covered_indices.insert(i);
        covered_indices.insert(j as usize);
    }

    result_edges
}

pub fn eulerian_tour(edges: &Vec<(usize, usize)>) -> Vec<usize> {
    /// Compute an Eulerian tour on the multigraph specified by the given edges.
    /// Edges are assumed to be undirected.
    ///
    /// It is assumed that every vertex has even degree.
    ///
    /// The result is returned as a sequence of vertices.
    let mut n_vertices = 0;
    for &edge in edges.iter() {
        if n_vertices <= edge.0 {
            n_vertices = edge.0 + 1;
        }
        if n_vertices <= edge.1 {
            n_vertices = edge.1 + 1;
        }
    }

    let mut neighborhoods: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
    for (i, &edge) in edges.iter().enumerate() {
        neighborhoods[edge.0].push(i);
        neighborhoods[edge.1].push(i);
    }

    let mut result = Vec::new();
    let mut stack = Vec::with_capacity(n_vertices);
    stack.push(edges[0].0);
    while stack.len() > 0 {
        let v = *stack.last().unwrap();
        if neighborhoods[v].len() == 0 {
            result.push(stack.pop().unwrap());
        } else {
            let edge_idx = neighborhoods[v].pop().unwrap();
            let edge = edges[edge_idx];
            let w = if edge.0 == v { edge.1 } else { edge.0 };
            let ew = neighborhoods[w]
                .iter()
                .position(|i| *i == edge_idx)
                .unwrap();
            neighborhoods[w].swap_remove(ew);
            stack.push(w);
        }
    }

    result
}

pub fn eulerian_tour_check(edges: &Vec<(usize, usize)>) -> Option<Vec<usize>> {
    /// Same as above, but we check that every vertex has even degree.

    let mut n_vertices = 0;
    for &edge in edges.iter() {
        if n_vertices <= edge.0 {
            n_vertices = edge.0 + 1;
        }
        if n_vertices <= edge.1 {
            n_vertices = edge.1 + 1;
        }
    }

    let mut neighborhoods: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
    for (i, &edge) in edges.iter().enumerate() {
        neighborhoods[edge.0].push(i);
        neighborhoods[edge.1].push(i);
    }
    for nbrhd in neighborhoods.iter() {
        //println!("{:?}", nbrhd.clone());
        if nbrhd.len() % 2 == 1 {
            return None;
        }
    }

    let mut result = Vec::new();
    let mut stack = Vec::with_capacity(n_vertices);
    stack.push(edges[0].0);
    while stack.len() > 0 {
        let v = *stack.last().unwrap();
        if neighborhoods[v].len() == 0 {
            result.push(stack.pop().unwrap());
        } else {
            let edge_idx = neighborhoods[v].pop().unwrap();
            let edge = edges[edge_idx];
            let w = if edge.0 == v { edge.1 } else { edge.0 };
            let ew = neighborhoods[w]
                .iter()
                .position(|i| *i == edge_idx)
                .unwrap();
            neighborhoods[w].swap_remove(ew);
            stack.push(w);
        }
    }

    Some(result)
}

fn eulerian_tour2<W: Copy>(
    edges: &Vec<(usize, usize, W)>,
    n_vertices: usize,
    neighborhoods: &Vec<Vec<usize>>,
) -> Vec<usize> {
    /// Same as `eulerian_tour` but the number of vertices is provided as an argument.
    let mut neighborhoods = neighborhoods.clone();
    let mut result = Vec::new();
    let mut stack = Vec::with_capacity(n_vertices);
    stack.push(edges[0].0);
    while stack.len() > 0 {
        let v = *stack.last().unwrap();
        if neighborhoods[v].len() == 0 {
            result.push(stack.pop().unwrap());
        } else {
            let edge_idx = neighborhoods[v].pop().unwrap();
            let edge = edges[edge_idx];
            let w = if edge.0 == v { edge.1 } else { edge.0 };
            let ew = neighborhoods[w]
                .iter()
                .position(|i| *i == edge_idx)
                .unwrap();
            neighborhoods[w].swap_remove(ew);
            stack.push(w);
        }
    }

    result
}

pub fn postman<W: Debug + Bounded + PartialOrd + Copy + Sized + Num + NumCast + Into<f64>>(
    weighted_edges: &Vec<(usize, usize, W)>,
) -> Vec<usize> {
    /// Solves the Chinese Postman Problem using a minimum weight T-join
    /// where T is the set of odd-degree vertices.
    ///
    /// Returns the postman tour as a sequence of vertices.

    let mut n_vertices = 0;
    for &edge in weighted_edges.iter() {
        if n_vertices <= edge.0 {
            n_vertices = edge.0 + 1;
        }
        if n_vertices <= edge.1 {
            n_vertices = edge.1 + 1;
        }
    }

    let mut neighborhoods: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
    for (i, &edge) in weighted_edges.iter().enumerate() {
        neighborhoods[edge.0].push(i);
        neighborhoods[edge.1].push(i);
    }
    let odd_verts: Vec<usize> = (0..n_vertices)
        .filter(|&v| neighborhoods[v].len() % 2 == 1)
        .collect();
    println!("odd_verts {:?}", odd_verts);

    if odd_verts.len() == 0 {
        eulerian_tour2(&weighted_edges, n_vertices, &neighborhoods)
    } else {
        let t_join = min_weight_t_join2(weighted_edges, &odd_verts, n_vertices);
        println!("t_join {:?}", Vec::from_iter(t_join.iter()));
        let new_edges = weighted_edges
            .iter()
            .map(|e| (e.0, e.1, ()))
            .chain(t_join.iter().map(|e| (e.0, e.1, ())))
            .collect();
        //println!("new_edges {:?}", new_edges);
        for (i, &new_edge) in t_join.iter().enumerate() {
            neighborhoods[new_edge.0].push(i + weighted_edges.len());
            neighborhoods[new_edge.1].push(i + weighted_edges.len());
        }
        eulerian_tour2(&new_edges, n_vertices, &neighborhoods)
    }
}
