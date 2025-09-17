use num::Bounded;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt::Display;

#[derive(PartialEq, Debug, Clone, Copy)]
enum UpdateMode {
    Undetermined,
    Vertex,
    SVertexFreeVertex,
    SBlossom,
    TBlossom,
}

pub struct BlossomData {
    n_vertices: usize,
    n_edges: usize,

    // label 0 corresponds to "unlabelled"
    // the search for augmenting paths / blossoms
    // begins at vertices with label 1
    // when an unlabelled vertex is reached from
    // a "label 1" vertex, label that vertex 2
    // likewise, when an unlabelled vertex is
    // reached from a "label 2" vertex, label
    // that vertex with 1
    pub blossom_labels: Vec<i8>,

    // if b is a top-level blossom
    // label_endpoints[b] records the vertex on
    // the other side of the edge through which b
    // received its label (i.e. the vertex that was
    // directly preceding blossom b in the search)
    pub label_endpoints: Vec<i64>,

    // if b is a non-trivial (non-vertex) sub-blossom
    blossom_endpoints: Vec<Vec<i64>>,

    // records the top-level blossom containing a given vertex
    pub blossom_id: Vec<usize>,

    // tracks the nesting of the at most
    // 2*n_vertices blossoms
    blossom_parent: Vec<i64>,
    blossom_children: Vec<Vec<usize>>,
    blossom_root: Vec<i64>,

    // for a free vertex or unreached
    // vertex inside a label 2 blossom,
    // this records the index of the edge
    // with least slack, or -1 if no such
    // edge exists
    pub best_edge: Vec<i64>,

    // for non-trivial (i.e. non-vertex) S-blossoms
    // this records the least-slack edges to neighboring
    // S-blossoms, which is useful for quickly computing
    // updates to the dual solution
    blossom_best_edges: Vec<Vec<usize>>,

    // track IDs which have not yet been used
    // to describe blossoms in the blossom tree
    unused_blossoms: Vec<usize>,

    // dual solution to the current matching
    // the matching is optimal when the slack is zero
    dual_soln: Vec<f64>,

    // allowed_edge[e] = true means that edge e
    // has zero slack, allowed_edge[e] = false
    // means the edge may or may not have zero slack
    pub allowed_edge: Vec<bool>,

    // records which way to
    update_mode: UpdateMode,
}

impl BlossomData {
    pub fn new(n_vertices: usize, n_edges: usize, max_weight: f64) -> Self {
        let blossom_labels = vec![0i8; 2 * n_vertices];
        let label_endpoints = vec![-1i64; 2 * n_vertices];
        let blossom_endpoints = vec![Vec::new(); 2 * n_vertices];
        let blossom_id = Vec::from_iter(0usize..n_vertices);
        let blossom_parent = vec![-1i64; 2 * n_vertices];
        let blossom_children = vec![Vec::new(); 2 * n_vertices];
        let mut blossom_root = Vec::from_iter(0i64..(n_vertices as i64));
        for _ in 0..n_vertices {
            blossom_root.push(-1);
        }
        let best_edge = vec![-1i64; 2 * n_vertices];
        let blossom_best_edges = vec![Vec::new(); 2 * n_vertices];
        let unused_blossoms = Vec::from_iter(n_vertices..2 * n_vertices);
        let mut dual_soln = vec![max_weight; n_vertices];
        for _ in 0..n_vertices {
            dual_soln.push(0.0);
        }
        let allowed_edge: Vec<bool> = vec![false; n_edges];
        let update_mode = UpdateMode::Undetermined;
        BlossomData {
            n_vertices,
            n_edges,
            blossom_labels,
            label_endpoints,
            blossom_endpoints,
            blossom_id,
            blossom_parent,
            blossom_children,
            blossom_root,
            best_edge,
            blossom_best_edges,
            unused_blossoms,
            dual_soln,
            allowed_edge,
            update_mode,
        }
    }

    pub fn clear(&mut self) {
        self.blossom_labels = vec![0; 2 * self.n_vertices];
        self.best_edge = vec![-1; 2 * self.n_vertices];
        self.blossom_best_edges = vec![Vec::new(); 2 * self.n_vertices];
        self.allowed_edge = vec![false; self.n_edges];
    }

    pub fn slack(&self, edge_idx: usize, weighted_edges: &Vec<(usize, usize, f64)>) -> f64 {
        let (i, j, w) = weighted_edges[edge_idx];
        self.dual_soln[i] + self.dual_soln[j] - 2.0 * w
    }

    fn compute_delta_vertices(&self) -> f64 {
        let mut d = f64::max_value();
        for i in 0..self.n_vertices {
            if self.dual_soln[i] < d {
                d = self.dual_soln[i];
            }
        }
        d
    }

    fn compute_delta_s_vertex_free_vertex(
        &mut self,
        delta: f64,
        weighted_edges: &Vec<(usize, usize, f64)>,
    ) -> (f64, usize) {
        let mut d = delta;
        let mut ix = 0usize;
        for v in 0..self.n_vertices {
            if self.blossom_labels[self.blossom_id[v]] == 0 && self.best_edge[v] != -1 {
                let de = self.slack(self.best_edge[v] as usize, weighted_edges);
                if de < d || self.update_mode == UpdateMode::Undetermined {
                    d = de;
                    ix = self.best_edge[v] as usize;
                    self.update_mode = UpdateMode::SVertexFreeVertex;
                }
            }
        }
        (d, ix)
    }

    fn compute_delta_s_blossoms(
        &mut self,
        delta: f64,
        best_edge_idx: usize,
        weighted_edges: &Vec<(usize, usize, f64)>,
    ) -> (f64, usize) {
        let mut d = delta;
        let mut ix = best_edge_idx;
        for b in 0..2 * self.n_vertices {
            if self.blossom_parent[b] == -1
                && self.blossom_labels[b] == 1
                && self.best_edge[b] != -1
            {
                let edge_idx = self.best_edge[b] as usize;
                let best_slack = 0.5 * self.slack(edge_idx, weighted_edges);
                if best_slack < d || self.update_mode == UpdateMode::Undetermined {
                    d = best_slack;
                    ix = edge_idx;
                    self.update_mode = UpdateMode::SBlossom;
                }
            }
        }
        (d, ix)
    }

    fn compute_delta_t_blossoms(&mut self, delta: f64) -> (usize, f64) {
        let mut d = delta;
        let mut blossom = 0;
        for b in self.n_vertices..2 * self.n_vertices {
            if self.blossom_root[b] >= 0
                && self.blossom_parent[b] == -1
                && self.blossom_labels[b] == 2
            {
                if self.dual_soln[b] < d || self.update_mode == UpdateMode::Undetermined {
                    d = self.dual_soln[b];
                    blossom = b;
                    self.update_mode = UpdateMode::TBlossom;
                }
            }
        }
        (blossom, d)
    }

    pub fn determine_delta_and_update_mode(
        &mut self,
        weighted_edges: &Vec<(usize, usize, f64)>,
        max_cardinality: bool,
    ) -> (f64, usize, usize) {
        let mut delta = <f64 as Bounded>::min_value();
        if !max_cardinality {
            delta = self.compute_delta_vertices();
            self.update_mode = UpdateMode::Vertex;
        }
        let (delta, best_edge) = self.compute_delta_s_vertex_free_vertex(delta, weighted_edges);
        let (delta, best_edge) = self.compute_delta_s_blossoms(delta, best_edge, weighted_edges);
        let (update_blossom, mut delta) = self.compute_delta_t_blossoms(delta);

        if self.update_mode == UpdateMode::Undetermined {
            //assert!(max_cardinality);
            // if we reach this point without determining an update mode
            // no more updates are necessary
            // in that case, mark the update mode as Vertex to terminate the loop after
            // one more update to the blossom structure
            self.update_mode = UpdateMode::Vertex;
            delta = self.compute_delta_vertices();
            delta = if delta < 0.0 { 0.0 } else { delta };
        }

        (delta, best_edge, update_blossom)
    }

    pub fn update_blossom_structure(
        &mut self,
        stack: &mut Vec<usize>,
        best_edge: usize,
        update_blossom: usize,
        weighted_edges: &Vec<(usize, usize, f64)>,
        endpoints: &Vec<usize>,
        matching: &Vec<i64>,
    ) -> bool {
        match self.update_mode {
            UpdateMode::Undetermined => {
                self.update_mode = UpdateMode::Undetermined;
            }
            UpdateMode::Vertex => {
                self.update_mode = UpdateMode::Undetermined;
                return true;
            }
            UpdateMode::SVertexFreeVertex => {
                self.update_mode = UpdateMode::Undetermined;
                self.allowed_edge[best_edge] = true;
                let (mut i, mut j, _) = weighted_edges[best_edge];
                if self.blossom_labels[self.blossom_id[i]] == 0 {
                    (j, i) = (i, j);
                }
                //assert_eq!(self.blossom_labels[self.blossom_id[i]], 1);
                stack.push(i);
            }
            UpdateMode::SBlossom => {
                self.update_mode = UpdateMode::Undetermined;
                self.allowed_edge[best_edge] = true;
                let (i, j, _) = weighted_edges[best_edge];
                //assert_eq!(self.blossom_labels[self.blossom_id[i]], 1);
                stack.push(i);
            }
            UpdateMode::TBlossom => {
                self.update_mode = UpdateMode::Undetermined;
                self.expand_blossom(stack, update_blossom, false, endpoints, matching);
            }
        }
        false
    }

    pub fn update_dual_soln(&mut self, delta: f64) {
        // the dual solution starts by satisfying the constraints
        // 1) For every i, j: dual_soln[i], dual_soln[j], slack(i, j) >= 0
        // 2) (i, j) is matched => slack(i, j) = 0
        // 3) dual_soln[b] > 0 => the maximum number of edges within blossom b have been matched
        // and we always update the dual solution so that it continues to
        // satsify (1) and (2)
        // the idea of the primal-dual algorithm is to gradually
        // reduce the number of violations to the 4th constraint
        // 4) i is single => dual_soln[i] = 0
        for v in 0..self.n_vertices {
            let label = self.blossom_labels[self.blossom_id[v]];
            if label == 1 {
                self.dual_soln[v] -= delta;
            } else if label == 2 {
                self.dual_soln[v] += delta;
            }
        }
        for b in self.n_vertices..2 * self.n_vertices {
            if self.blossom_root[b] >= 0 && self.blossom_parent[b] == -1 {
                if self.blossom_labels[b] == 1 {
                    self.dual_soln[b] += delta;
                } else if self.blossom_labels[b] == 2 {
                    self.dual_soln[b] -= delta;
                }
            }
        }
    }

    pub fn s_blossom_is_tight(&mut self, blossom_id: usize) -> bool {
        self.blossom_parent[blossom_id] == -1
            && self.blossom_root[blossom_id] >= 0
            && self.blossom_labels[blossom_id] == 1
            && self.dual_soln[blossom_id] < 1e-12
    }

    pub fn collect_leaves(&self, current_blossom: usize, deep: bool) -> Vec<usize> {
        if current_blossom < self.n_vertices {
            return vec![current_blossom];
        }

        // another heuristic to balance cache hits vs heap allocations
        let alloc = if !deep {
            2 * self.blossom_children[current_blossom].len()
        } else {
            3
        };
        let mut leaves = Vec::with_capacity(alloc);

        for &b in self.blossom_children[current_blossom].iter() {
            if b < self.n_vertices {
                leaves.push(b)
            } else {
                leaves.append(&mut self.collect_leaves(b, true));
            }
        }
        leaves
    }

    pub fn assign_label<'a>(
        &'a mut self,
        stack: &mut Vec<usize>,
        vertex: usize,
        label: u8,
        endpoint: i64,
        endpoints: &Vec<usize>,
        matching: &Vec<i64>,
    ) {
        //println!("assign_label: {} {} {}", vertex, label, endpoint);
        let b = self.blossom_id[vertex];
        //assert_eq!(self.blossom_labels[vertex], 0);
        //assert_eq!(self.blossom_labels[b], 0);
        self.blossom_labels[vertex] = label as i8;
        self.blossom_labels[b] = label as i8;
        self.label_endpoints[vertex] = endpoint;
        self.label_endpoints[b] = endpoint;
        self.best_edge[vertex] = -1;
        self.best_edge[b] = -1;
        if label == 1 {
            ////println!("Append: {:?}", self.collect_leaves(b, false));
            stack.append(&mut self.collect_leaves(b, false));
        } else if label == 2 {
            let root = self.blossom_root[b] as usize;
            //assert!(matching[root] >= 0);
            ////println!("label endpoints {:?}", self.label_endpoints);
            self.assign_label(
                stack,
                endpoints[matching[root] as usize],
                1,
                matching[root] ^ 1,
                endpoints,
                matching,
            )
        } else {
            panic!("Tried to assign label other than 1 or 2");
        }
    }

    pub fn check_for_blossom_or_augmenting_path(
        &mut self,
        v: usize,
        w: usize,
        endpoints: &Vec<usize>,
        matching: &Vec<i64>,
    ) -> Option<usize> {
        //println!("SCAN {} {}", v, w);
        let mut v = v as i64;
        let mut w = w as i64;
        let mut path = Vec::new();
        while v != -1 || w != -1 {
            let mut b = self.blossom_id[v as usize];
            // we label traversed vertices by 3,
            // so that we can detect an odd cycle (blossom)
            // when we encounter them again we know we have
            // detected an odd cycle (blossom)
            if self.blossom_labels[b] > 2 {
                for c in path.iter() {
                    self.blossom_labels[*c] = 1;
                }
                return Some(self.blossom_root[b] as usize);
            }

            //assert_eq!(self.blossom_labels[b], 1);
            path.push(b);
            // label 3 here so we know which
            // vertices we have already traced
            self.blossom_labels[b] = 3;
            //assert_eq!(
            //    self.label_endpoints[b],
            //    matching[self.blossom_root[b] as usize]
            //);
            // if no blossom is detected, trace our
            // path through the search forest backwards
            if self.label_endpoints[b] == -1 {
                v = -1;
            } else {
                v = endpoints[self.label_endpoints[b] as usize] as i64;
                b = self.blossom_id[v as usize];
                //assert_eq!(self.blossom_labels[b], 2);
                //assert!(self.label_endpoints[b] >= 0);
                v = endpoints[self.label_endpoints[b] as usize] as i64;
            }
            // and alternate the paths which we trace back
            if w != -1 {
                (v, w) = (w, v);
            }
        }
        for &b in path.iter() {
            self.blossom_labels[b] = 1;
        }
        None
    }

    pub fn add_blossom<'a>(
        &'a mut self,
        stack: &mut Vec<usize>,
        root: usize,
        edge_idx: usize,
        matching: &Vec<i64>,
        endpoints: &Vec<usize>,
        nbrhd_endpoints: &Vec<Vec<usize>>,
        weighted_edges: &Vec<(usize, usize, f64)>,
    ) {
        let (mut v, mut w, _) = weighted_edges[edge_idx];
        let bb = self.blossom_id[root];
        let mut bv = self.blossom_id[v];
        let mut bw = self.blossom_id[w];
        let b = match self.unused_blossoms.pop() {
            None => panic!("No more unused blossoms."),
            Some(c) => c,
        };
        //println!("ADD BLOSSOM {} {} {} {} {}", root, edge_idx, v, w, b);

        self.blossom_root[b] = root as i64;
        self.blossom_parent[b] = -1;
        self.blossom_parent[bb] = b as i64;

        self.blossom_children[b] = Vec::new();
        self.blossom_endpoints[b] = Vec::new();

        while bv != bb {
            self.blossom_parent[bv] = b as i64;
            self.blossom_children[b].push(bv);
            self.blossom_endpoints[b].push(self.label_endpoints[bv]);

            v = endpoints[self.label_endpoints[bv] as usize];
            bv = self.blossom_id[v];
        }

        self.blossom_children[b].push(bb);
        self.blossom_children[b] = self.blossom_children[b].iter().copied().rev().collect();
        self.blossom_endpoints[b] = self.blossom_endpoints[b].iter().copied().rev().collect();
        self.blossom_endpoints[b].push(2 * edge_idx as i64);
        //println!("BLOSSOM ENDPOINT {} PUSH {}", b, 2*edge_idx as i64);

        while bw != bb {
            self.blossom_parent[bw] = b as i64;
            self.blossom_children[b].push(bw);
            self.blossom_endpoints[b].push(self.label_endpoints[bw] ^ 1);

            w = endpoints[self.label_endpoints[bw] as usize];
            bw = self.blossom_id[w];
        }

        self.blossom_labels[b] = 1;
        self.label_endpoints[b] = self.label_endpoints[bb];
        self.dual_soln[b] = 0.0;
        let leaves: Vec<usize> = self.collect_leaves(b, false);
        let mut two_leaves = Vec::with_capacity(leaves.len());
        for &leaf in leaves.iter() {
            if self.blossom_labels[self.blossom_id[leaf]] == 2 {
                two_leaves.push(leaf);
            }
            self.blossom_id[leaf] = b;
        }

        let mut best_edges = vec![-1i64; 2 * self.n_vertices];
        for &bv in self.blossom_children[b].iter() {
            if self.blossom_best_edges[bv].len() == 0 {
                // if the blossom has not stored its least-slack edges, iterate through all of its vertices
                for nbr_list in self
                    .collect_leaves(bv, false)
                    .iter()
                    .map(|v| nbrhd_endpoints[*v].iter().map(|p| *p / 2))
                {
                    for nbr in nbr_list {
                        let (mut i, mut j, _) = weighted_edges[nbr];
                        if self.blossom_id[j] == b {
                            (j, i) = (i, j);
                        }
                        let bj = self.blossom_id[j];
                        if bj != b
                            && self.blossom_labels[bj] == 1
                            && (best_edges[bj] == -1
                                || self.slack(nbr, weighted_edges)
                                    < self.slack(best_edges[bj] as usize, weighted_edges))
                        {
                            best_edges[bj] = nbr as i64;
                        }
                    }
                }
            } else {
                // otherwise we can iterate through the blossom's pre-computed least slack edges
                for nbr in self.blossom_best_edges[bv].iter() {
                    let (mut i, mut j, _) = weighted_edges[*nbr];
                    if self.blossom_id[j] == b {
                        (j, i) = (i, j);
                    }
                    let bj = self.blossom_id[j];
                    if bj != b
                        && self.blossom_labels[bj] == 1
                        && (best_edges[bj] == -1
                            || self.slack(*nbr, weighted_edges)
                                < self.slack(best_edges[bj] as usize, weighted_edges))
                    {
                        best_edges[bj] = *nbr as i64;
                    }
                }
            }
            self.blossom_best_edges[bv] = Vec::new();
            self.best_edge[bv] = -1;
        }
        self.blossom_best_edges[b] = best_edges
            .iter()
            .filter(|k| **k > -1)
            .map(|k| *k as usize)
            .collect();
        self.best_edge[b] = -1;
        for k in self.blossom_best_edges[b].iter() {
            if self.best_edge[b] == -1
                || self.slack(*k, weighted_edges)
                    < self.slack(self.best_edge[b] as usize, weighted_edges)
            {
                self.best_edge[b] = *k as i64;
            }
        }

        //println!("APPENDING {:?}", two_leaves);
        stack.append(&mut two_leaves);
        //println!("BLOSSOM {:?}", self.blossom_children);
    }

    pub fn expand_blossom<'a>(
        &'a mut self,
        stack: &mut Vec<usize>,
        blossom_id: usize,
        endstage: bool,
        endpoints: &Vec<usize>,
        matching: &Vec<i64>,
    ) {
        let children: Vec<usize> = self.blossom_children[blossom_id].clone();
        for &c in children.iter() {
            self.blossom_parent[c] = -1;
            if c < self.n_vertices {
                self.blossom_id[c] = c
            } else if endstage && self.dual_soln[c] < 1e-12 {
                self.expand_blossom(stack, c, endstage, endpoints, matching);
            } else {
                let leaves: Vec<usize> = self.collect_leaves(c, false);
                for v in leaves {
                    self.blossom_id[v] = c;
                }
            }
        }

        if !endstage && self.blossom_labels[blossom_id] == 2 {
            //assert!(self.label_endpoints[blossom_id] >= 0);
            let entry_child =
                self.blossom_id[endpoints[(self.label_endpoints[blossom_id] ^ 1) as usize]];
            let j = match self.blossom_children[blossom_id]
                .iter()
                .position(|v| *v == entry_child)
            {
                None => panic!("Did not find entry_child in children"),
                Some(x) => x,
            };

            let mut endpoint = self.label_endpoints[blossom_id];
            if j % 2 == 1 {
                let mut k = j;
                while k < self.blossom_children[blossom_id].len() {
                    let ix = (endpoint ^ 1) as usize;
                    self.blossom_labels[endpoints[ix]] = 0;
                    self.blossom_labels
                        [endpoints[(self.blossom_endpoints[blossom_id][k] ^ 1) as usize]] = 0;
                    self.assign_label(stack, endpoints[ix], 2, endpoint, endpoints, matching);
                    self.allowed_edge[(self.blossom_endpoints[blossom_id][k] / 2) as usize] = true;
                    k += 1;
                    endpoint = self.blossom_endpoints[blossom_id][k];
                    self.allowed_edge[(endpoint / 2) as usize] = true;
                    k += 1;
                }

                k = 0;
                let mut bv = self.blossom_children[blossom_id][k];
                let ix = (endpoint ^ 1) as usize;
                self.blossom_labels[endpoints[ix]] = 2;
                self.blossom_labels[bv] = 2;
                self.label_endpoints[endpoints[ix]] = endpoint;
                self.label_endpoints[bv] = endpoint;
                self.best_edge[bv] = -1;

                k = 1;
                while self.blossom_children[blossom_id][j] != entry_child {
                    bv = self.blossom_children[blossom_id][j];
                    if self.blossom_labels[bv] == 1 {
                        k += 1;
                    } else {
                        let mut v = None;
                        for &leaf in self.collect_leaves(bv, false).iter() {
                            if self.blossom_labels[leaf] != 0 {
                                v = Some(leaf);
                                break;
                            }
                        }
                        match v {
                            None => continue,
                            Some(v) => {
                                //assert_eq!(self.blossom_labels[v], 2);
                                //assert_eq!(self.blossom_id[v], bv);
                                self.blossom_labels[v] = 0;
                                self.blossom_labels[endpoints
                                    [matching[self.blossom_root[bv] as usize] as usize]] = 0;
                                self.assign_label(stack, v, 2, endpoint, endpoints, matching);
                            }
                        }
                        k += 1;
                    }
                }
            } else {
                let mut k = j;
                while k > 0 {
                    let ix = (endpoint ^ 1) as usize;
                    self.blossom_labels[endpoints[ix]] = 0;
                    self.blossom_labels
                        [endpoints[(self.blossom_endpoints[blossom_id][k - 1]) as usize]] = 0;
                    self.assign_label(stack, endpoints[ix], 2, endpoint, endpoints, matching);
                    self.allowed_edge[(self.blossom_endpoints[blossom_id][k - 1] / 2) as usize] =
                        true;
                    k -= 1;
                    endpoint = self.blossom_endpoints[blossom_id][k - 1] ^ 1;
                    self.allowed_edge[(endpoint / 2) as usize] = true;
                    k -= 1;
                }

                let mut bv = self.blossom_children[blossom_id][0];
                let ix = (endpoint ^ 1) as usize;
                self.blossom_labels[endpoints[ix]] = 2;
                self.blossom_labels[bv] = 2;
                self.label_endpoints[endpoints[ix]] = endpoint;
                self.label_endpoints[bv] = endpoint;
                self.best_edge[bv] = -1;

                k = self.blossom_children[blossom_id].len() - 1;
                while self.blossom_children[blossom_id][j] != entry_child {
                    bv = self.blossom_children[blossom_id][j];
                    if self.blossom_labels[bv] == 1 {
                        k -= 1;
                    } else {
                        let mut v = None;
                        for &leaf in self.collect_leaves(bv, false).iter() {
                            if self.blossom_labels[leaf] != 0 {
                                v = Some(leaf);
                                break;
                            }
                        }
                        match v {
                            None => continue,
                            Some(v) => {
                                //assert_eq!(self.blossom_labels[v], 2);
                                //assert_eq!(self.blossom_id[v], bv);
                                self.blossom_labels[v] = 0;
                                self.blossom_labels[endpoints
                                    [matching[self.blossom_root[bv] as usize] as usize]] = 0;
                                self.assign_label(stack, v, 2, endpoint, endpoints, matching);
                            }
                        }
                        k -= 1;
                    }
                }
            }
        }

        // recycle the blossom id
        self.blossom_labels[blossom_id] = -1;
        self.label_endpoints[blossom_id] = -1;
        self.blossom_root[blossom_id] = -1;
        self.blossom_best_edges[blossom_id] = Vec::new();
        self.best_edge[blossom_id] = -1;
        self.unused_blossoms.push(blossom_id);
    }

    // this function augments the matching within and around the blossom
    pub fn augment_blossom(
        &mut self,
        blossom_id: usize,
        vertex: usize,
        endpoints: &Vec<usize>,
        matching: &mut Vec<i64>,
    ) {
        ////println!("{:?}", self.blossom_endpoints[blossom_id]);
        // let t be the blossom which contains vertex
        // and which is a direct child of blossom_id
        let mut t = vertex;
        while self.blossom_parent[t] as usize != blossom_id {
            ////println!("{}", t);
            t = self.blossom_parent[t] as usize;
        }
        // if t is actually a non-trivial (non-vertex) blossom
        // we expand it before continuing
        if t >= self.n_vertices {
            self.augment_blossom(t, vertex, endpoints, matching);
        }

        // now find the index of t within the children of blossom_id
        let i = self.blossom_children[blossom_id]
            .iter()
            .position(|c| *c == t)
            .unwrap();
        //println!("AUGMENT BLOSSOM {blossom_id} {vertex} {i}");
        let mut j = i;
        if i % 2 == 1 {
            ////println!("{}", self.blossom_children[blossom_id].len());
            while j != 0 {
                j += 1;
                //println!("{}", j);
                t = self.blossom_children[blossom_id][j];
                let p = self.blossom_endpoints[blossom_id][j] as usize;
                let p_parity = p ^ 1;
                ////println!("{} {} {}", j, p, i%2);
                if t >= self.n_vertices {
                    ////println!("1142");
                    self.augment_blossom(t, endpoints[p], endpoints, matching);
                }

                j += 1;
                j = j % self.blossom_children[blossom_id].len();
                //println!("{}", j);
                ////println!("{} {} {}", j, p, p_parity);
                t = self.blossom_children[blossom_id][j];
                if t >= self.n_vertices {
                    ////println!("1151");
                    ////println!("{:?}", endpoints);
                    self.augment_blossom(t, endpoints[p_parity], endpoints, matching);
                }
                //println!("BLOSSOM MATCHING ASSIGNMENT {} {} {} {}", endpoints[p], p_parity, endpoints[p_parity], p);
                matching[endpoints[p]] = p_parity as i64;
                matching[endpoints[p_parity]] = p as i64;
            }
        } else {
            let last_ix = self.blossom_children[blossom_id].len() - 1;
            while j != 0 {
                j -= 1;
                //println!("{}", j);
                t = self.blossom_children[blossom_id][j];
                let p = (self.blossom_endpoints[blossom_id][j - 1] ^ 1) as usize;
                let p_parity = p ^ 1;
                if t >= self.n_vertices {
                    ////println!("1166");
                    self.augment_blossom(t, endpoints[p], endpoints, matching);
                }

                j = if j == 0 { last_ix } else { j - 1 };
                //println!("{}", j);
                t = self.blossom_children[blossom_id][j];
                if t >= self.n_vertices {
                    ////println!("1173");
                    self.augment_blossom(t, endpoints[p_parity], endpoints, matching);
                }

                //println!("BLOSSOM MATCHING ASSIGNMENT {} {} {} {}", endpoints[p], p_parity, endpoints[p_parity], p);
                matching[endpoints[p]] = p_parity as i64;
                matching[endpoints[p_parity]] = p as i64;
            }
        }

        self.blossom_children[blossom_id].rotate_left(i);
        self.blossom_endpoints[blossom_id].rotate_left(i);
        self.blossom_root[blossom_id] = self.blossom_root[self.blossom_children[blossom_id][0]];
        //assert_eq!(self.blossom_root[blossom_id] as usize, vertex);
    }

    pub fn augment_matching(
        &mut self,
        edge_idx: usize,
        weighted_edges: &Vec<(usize, usize, f64)>,
        endpoints: &Vec<usize>,
        matching: &mut Vec<i64>,
    ) {
        let (v, w, _) = weighted_edges[edge_idx];
        let (mut s, mut p) = (v, 2 * edge_idx + 1);
        ////println!("AUGMENT MATCHING {} {} {}", edge_idx, v, w);
        loop {
            let bs = self.blossom_id[s];
            ////println!("blossom_labels {:?}", self.blossom_labels.clone());
            //assert_eq!(self.blossom_labels[bs], 1);
            //assert_eq!(
            //    self.label_endpoints[bs],
            //    matching[self.blossom_root[bs] as usize]
            //);
            if bs >= self.n_vertices {
                self.augment_blossom(bs, s, endpoints, matching);
            }
            //println!("A MATE {} ASSIGNED {}", s, p);
            matching[s] = p as i64;

            if self.label_endpoints[bs] == -1 {
                break;
            }

            let t = endpoints[self.label_endpoints[bs] as usize];
            let bt = self.blossom_id[t];
            //assert_eq!(self.blossom_labels[bt], 2);
            //assert!(self.label_endpoints[bt] >= 0);
            s = endpoints[self.label_endpoints[bt] as usize];
            let j = endpoints[(self.label_endpoints[bt] as usize) ^ 1];
            //assert_eq!(self.blossom_root[bt] as usize, t);
            if bt >= self.n_vertices {
                self.augment_blossom(bt, j, endpoints, matching);
            }
            //println!("B MATE {} ASSIGNED {}", j, self.label_endpoints[bt]);
            matching[j] = self.label_endpoints[bt];
            p = (self.label_endpoints[bt] ^ 1) as usize;
        }

        let (mut s, mut p) = (w, 2 * edge_idx);
        loop {
            let bs = self.blossom_id[s];
            //assert_eq!(self.blossom_labels[bs], 1);
            //assert_eq!(
            //    self.label_endpoints[bs],
            //    matching[self.blossom_root[bs] as usize]
            //);
            if bs >= self.n_vertices {
                self.augment_blossom(bs, s, endpoints, matching);
            }
            matching[s] = p as i64;
            //println!("C MATE {} ASSIGNED {}", s, p);

            if self.label_endpoints[bs] == -1 {
                break;
            }

            let t = endpoints[self.label_endpoints[bs] as usize];
            let bt = self.blossom_id[t];
            //assert_eq!(self.blossom_labels[bt], 2);
            //assert!(self.label_endpoints[bt] >= 0);
            s = endpoints[self.label_endpoints[bt] as usize];
            let j = endpoints[(self.label_endpoints[bt] as usize) ^ 1];
            //assert_eq!(self.blossom_root[bt] as usize, t);
            if bt >= self.n_vertices {
                self.augment_blossom(bt, j, endpoints, matching);
            }
            matching[j] = self.label_endpoints[bt];
            //println!("D MATE {} ASSIGNED {}", j, self.label_endpoints[bt]);
            p = (self.label_endpoints[bt] ^ 1) as usize;
        }
    }
}
