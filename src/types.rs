use num::{Bounded, One};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt::Display;
use std::iter::{Chain, Copied};
use std::slice::Iter;
use std::vec::IntoIter;

pub type Edge<W> = (usize, usize, W);

fn get_vertex_set(edges: &Vec<(usize, usize)>) -> BTreeSet<usize> {
    let mut vert_set = BTreeSet::new();
    for edge in edges.iter() {
        vert_set.insert(edge.0);
        vert_set.insert(edge.1);
    }
    vert_set
}

pub trait Neighborhood {
    fn contains(&self, item: usize) -> bool;
    fn insert(&mut self, item: usize);
    fn remove(&mut self, item: usize) -> bool;
    fn into_iter_no_move(&self) -> impl Iterator<Item = usize>;
    fn new() -> Self;
}

impl Neighborhood for Vec<usize> {
    fn contains(&self, item: usize) -> bool {
        self.iter().any(|x| *x == item)
    }

    fn insert(&mut self, item: usize) {
        if !self.contains(item) {
            self.push(item);
        }
    }

    fn remove(&mut self, item: usize) -> bool {
        match self.iter().position(|x| *x == item) {
            Some(i) => {
                self.remove(i);
                true
            }
            None => false,
        }
    }

    fn into_iter_no_move(&self) -> impl Iterator<Item = usize> {
        self.iter().map(|x| *x)
    }

    fn new() -> Self {
        Vec::new()
    }
}

impl Neighborhood for HashSet<usize> {
    fn contains(&self, item: usize) -> bool {
        self.contains(&item)
    }

    fn insert(&mut self, item: usize) {
        self.insert(item);
    }

    fn remove(&mut self, item: usize) -> bool {
        self.remove(&item)
    }

    fn into_iter_no_move(&self) -> impl Iterator<Item = usize> {
        self.iter().map(|x| *x)
    }

    fn new() -> Self {
        HashSet::new()
    }
}

impl Neighborhood for BTreeSet<usize> {
    fn contains(&self, item: usize) -> bool {
        self.contains(&item)
    }

    fn insert(&mut self, item: usize) {
        self.insert(item);
    }

    fn remove(&mut self, item: usize) -> bool {
        self.remove(&item)
    }

    fn into_iter_no_move(&self) -> impl Iterator<Item = usize> {
        self.iter().map(|x| *x)
    }

    fn new() -> Self {
        BTreeSet::new()
    }
}

pub struct SquareMatrix<T> {
    pub data: Vec<T>,
    pub len: usize,
}

impl<T: Sized + Copy> SquareMatrix<T> {
    pub fn from_weighted_edges(
        weighted_edges: &Vec<Edge<T>>,
        directed: bool,
        no_edge_value: T,
        num_vertices: Option<usize>,
    ) -> Self {
        if weighted_edges.len() == 0 {
            return SquareMatrix {
                data: Vec::new(),
                len: 0,
            };
        }

        let num_vertices = match num_vertices {
            None => get_vertex_set(&weighted_edges.iter().map(|e| (e.0, e.1)).collect()).len(),
            Some(n) => n,
        };

        let mut result_data = vec![no_edge_value; num_vertices * num_vertices];
        for (i, j, w) in weighted_edges.iter() {
            result_data[num_vertices * i + j] = *w;
            if !directed {
                result_data[num_vertices * j + i] = *w;
            }
        }

        SquareMatrix {
            data: result_data,
            len: num_vertices,
        }
    }

    pub fn get_ix(&self, i: usize, j: usize) -> T {
        self.data[self.len * i + j]
    }

    pub fn set_ix(&mut self, i: usize, j: usize, val: T) {
        self.data[self.len * i + j] = val;
    }
}

impl<T: Display> SquareMatrix<T> {
    fn intersperse_spaces(st: &[String]) -> String {
        let mut result = String::new();
        for s in st.iter().cloned() {
            result.push_str(&s);
            result.push(' ');
        }
        result
    }

    pub fn to_string(&self) -> String {
        let helper: Vec<String> = self
            .data
            .iter()
            .map(|x| format!("{}", x))
            .collect::<Vec<String>>()
            .chunks(self.len)
            .map(|s| Self::intersperse_spaces(s))
            .collect();
        let mut result = String::new();
        for s in helper.iter() {
            result.push_str(s);
            result.push('\n');
        }
        result
    }
}

pub struct Graph<W, N: Neighborhood + FromIterator<usize> + Clone> {
    // Set of vertices of the graph
    pub vertices: BTreeSet<usize>,
    // Set of adjacencies of the graph using some arbitrary neighborhood datastructure.
    // It is assumed that this vector is long enough such that the vertices may be used
    // as indices into the vector, all though not all vertices in the range 0..max vertex
    // may be present in the vertex set. Adjacencies corresponding to non-existent vertices
    // should be kept empty
    pub adjacency_list: Vec<N>,
    pub edges: Option<Vec<(usize, usize)>>,
    pub weight_matrix: Option<SquareMatrix<W>>,
    pub directed: bool,
}

impl<W, N: Neighborhood + FromIterator<usize> + Clone> Graph<W, N> {
    pub fn compute_edges(&mut self) {
        match &self.edges {
            None => {
                let mut es = BTreeSet::new();
                for v in self.vertices.iter() {
                    for w in self.adjacency_list[*v].into_iter_no_move() {
                        if self.directed || !es.contains(&(w, *v)) {
                            es.insert((*v, w));
                        }
                    }
                }
                self.edges = Some(Vec::from_iter(es.iter().map(|e| *e)));
            }
            Some(_) => return,
        }
    }

    pub fn get_edges_set(&self) -> BTreeSet<(usize, usize)> {
        match &self.edges {
            None => {
                let mut es = BTreeSet::new();
                for v in self.vertices.iter() {
                    for w in self.adjacency_list[*v].into_iter_no_move() {
                        if self.directed || !es.contains(&(w, *v)) {
                            es.insert((*v, w));
                        }
                    }
                }
                es
            }
            Some(es) => BTreeSet::from_iter(es.iter().map(|e| *e)),
        }
    }

    pub fn get_edges_vec(&self) -> Vec<(usize, usize)> {
        match &self.edges {
            None => {
                let mut es = BTreeSet::new();
                for v in self.vertices.iter() {
                    for w in self.adjacency_list[*v].into_iter_no_move() {
                        if self.directed || !es.contains(&(w, *v)) {
                            es.insert((*v, w));
                        }
                    }
                }
                Vec::from_iter(es.into_iter())
            }
            Some(es) => es.clone(),
        }
    }
}

impl<N: Neighborhood + FromIterator<usize> + Clone> Graph<(), N> {
    pub fn from_edges(edges: &Vec<(usize, usize)>, directed: bool) -> Self {
        let vertices = get_vertex_set(&edges);
        let mut adjacency_list = vec![N::new(); vertices.len()];
        if directed {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
            }
        } else {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
                adjacency_list[edge.1].insert(edge.0);
            }
        }
        //dbg!(adjacency_list.iter().map(|n| Vec::from_iter(n.into_iter_no_move())).collect::<Vec<Vec<usize>>>());
        let edges = Some(edges.iter().map(|e| (e.0, e.1)).collect());
        let weight_matrix = None;

        Graph {
            vertices,
            adjacency_list,
            edges,
            weight_matrix,
            directed,
        }
    }

    pub fn from_edges_and_vertices(
        edges: &Vec<(usize, usize)>,
        vertices: &Vec<usize>,
        directed: bool,
    ) -> Self {
        let vertices = BTreeSet::from_iter(vertices.iter().map(|v| *v));
        let edge_vertices = get_vertex_set(&edges);
        if !vertices.is_superset(&edge_vertices) {
            panic!("Passed bad vertex set!");
        }

        let mut adjacency_list = vec![N::new(); vertices.len()];
        if directed {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
            }
        } else {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
                adjacency_list[edge.1].insert(edge.0);
            }
        }
        //dbg!(adjacency_list.iter().map(|n| Vec::from_iter(n.into_iter_no_move())).collect::<Vec<Vec<usize>>>());
        let edges = Some(edges.iter().map(|e| (e.0, e.1)).collect());
        let weight_matrix = None;

        Graph {
            vertices,
            adjacency_list,
            edges,
            weight_matrix,
            directed,
        }
    }
}

impl<W: Bounded + Sized + Copy, N: Neighborhood + FromIterator<usize> + Clone> Graph<W, N> {
    pub fn from_weighted_edges(edges: &Vec<Edge<W>>, directed: bool) -> Self {
        let vertices = get_vertex_set(&edges.iter().map(|e| (e.0, e.1)).collect());
        let mut adjacency_list = vec![N::new(); vertices.len()];
        if directed {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
            }
        } else {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
                adjacency_list[edge.1].insert(edge.0);
            }
        }

        let weight_matrix = Some(SquareMatrix::from_weighted_edges(
            &edges,
            directed,
            W::max_value(),
            Some(vertices.len()),
        ));
        let edges = Some(edges.into_iter().map(|e| (e.0, e.1)).collect());

        Graph {
            vertices,
            adjacency_list,
            edges,
            weight_matrix,
            directed,
        }
    }

    pub fn from_weighted_edges_and_vertices(
        edges: &Vec<Edge<W>>,
        vertices: &Vec<usize>,
        directed: bool,
    ) -> Self {
        let vertices = BTreeSet::from_iter(vertices.iter().map(|v| *v));
        let edge_vertices = get_vertex_set(&edges.iter().map(|e| (e.0, e.1)).collect());
        if !vertices.is_superset(&edge_vertices) {
            panic!("Passed bad vertex set!");
        }

        let mut adjacency_list = vec![N::new(); vertices.len()];
        if directed {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
            }
        } else {
            for edge in edges.iter() {
                adjacency_list[edge.0].insert(edge.1);
                adjacency_list[edge.1].insert(edge.0);
            }
        }

        let weight_matrix = Some(SquareMatrix::from_weighted_edges(
            &edges,
            directed,
            W::max_value(),
            Some(vertices.len()),
        ));
        let edges = Some(edges.into_iter().map(|e| (e.0, e.1)).collect());

        Graph {
            vertices,
            adjacency_list,
            edges,
            weight_matrix,
            directed,
        }
    }
}

impl<W, N: Neighborhood + FromIterator<usize> + Clone> Graph<W, N> {
    pub fn check_for_unmarked_edge(
        &self,
        marked_edges: &HashSet<(usize, usize)>,
        incident: usize,
    ) -> Option<usize> {
        // find an unmarked edge which is incident to the given vertex, if one exists
        // return the vertex opposite to the given vertex
        for v in self.adjacency_list[incident].into_iter_no_move() {
            if !marked_edges.contains(&(incident, v)) {
                return Some(v);
            }
        }
        None
    }
}

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
    pub blossom_labels: Vec<u8>,

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

pub struct BlossomLeaves<'a> {
    blossom_data: &'a BlossomData,
    current_blossom: usize,
    branch_path: Vec<usize>,
    done: bool,
}

impl<'a> Iterator for BlossomLeaves<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            None
        } else if self.current_blossom < self.blossom_data.n_vertices {
            // climb back up the tree until we find a branch to our right which we did not yet descend
            while let Some(i) = self.branch_path.pop() {
                let temp = self.current_blossom;
                let parent = self.blossom_data.blossom_parent[temp];
                if parent < 0 {
                    // this condition should actually be redundant
                    // with branch_path.pop() being None
                    self.done = true;
                    return Some(temp);
                } else if i < self.blossom_data.blossom_children[parent as usize].len() - 1 {
                    // start descending this branch
                    self.current_blossom =
                        self.blossom_data.blossom_children[parent as usize][i + 1];
                    self.branch_path.push(i + 1);
                    return Some(temp);
                }
                // otherwise continue ascending the tree
                self.current_blossom = parent as usize;
            }
            self.done = true;
            Some(self.current_blossom)
        } else {
            self.branch_path.push(0);
            self.current_blossom = self.blossom_data.blossom_children[self.current_blossom][0];
            self.next()
        }
    }
}

impl<'a> BlossomLeaves<'a> {
    fn new(current_blossom: usize, blossom_data: &'a BlossomData) -> Self {
        // heuristic for minimizing heap allocations while maximizing cache hits
        // letting branch_path be the largest capacity of n_vertices would totally
        // minimize heap allocations, but it could be too big to fit in the cache
        let branch_path = Vec::with_capacity(8);
        BlossomLeaves {
            blossom_data: blossom_data,
            current_blossom: current_blossom,
            branch_path: branch_path,
            done: false,
        }
    }
}

impl BlossomData {
    pub fn new(n_vertices: usize, n_edges: usize) -> Self {
        let blossom_labels = vec![0u8; 2 * n_vertices];
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
        let mut dual_soln = vec![<f64 as Bounded>::max_value(); n_vertices];
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
        self.blossom_best_edges = vec![Vec::new(); self.n_vertices];
        self.allowed_edge = vec![false; self.n_edges];
        self.update_mode = UpdateMode::Undetermined;
    }

    pub fn slack(
        &self,
        edge_idx: usize,
        edges: &Vec<(usize, usize)>,
        weight_matrix: &SquareMatrix<f64>,
    ) -> f64 {
        let (i, j) = edges[edge_idx];
        self.dual_soln[i] + self.dual_soln[j] - 2.0 * weight_matrix.get_ix(i, j)
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
        edges: &Vec<(usize, usize)>,
        weight_matrix: &SquareMatrix<f64>,
    ) -> (f64, usize) {
        let mut d = delta;
        let mut ix = 0usize;
        for v in 0..self.n_vertices {
            if self.blossom_labels[self.blossom_id[v]] == 0 && self.best_edge[v] != -1 {
                let de = self.slack(self.best_edge[v] as usize, edges, weight_matrix);
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
        edges: &Vec<(usize, usize)>,
        weight_matrix: &SquareMatrix<f64>,
    ) -> (f64, usize) {
        let mut d = delta;
        let mut ix = best_edge_idx;
        for b in 0..2 * self.n_vertices {
            if self.blossom_parent[b] == -1
                && self.blossom_labels[b] == 1
                && self.best_edge[b] != -1
            {
                let edge_idx = self.best_edge[b] as usize;
                let best_slack = 0.5 * self.slack(edge_idx, edges, weight_matrix);
                if best_slack < d || self.update_mode == UpdateMode::Undetermined {
                    d = best_slack;
                    ix = edge_idx;
                    self.update_mode = UpdateMode::SBlossom;
                }
            }
        }
        (d, ix)
    }

    fn compute_delta_t_blossoms(
        &mut self,
        delta: f64
    ) -> (usize, f64) {
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
        edges: &Vec<(usize, usize)>,
        weight_matrix: &SquareMatrix<f64>,
        max_cardinality: bool,
    ) -> (f64, usize, usize) {
        let mut delta = <f64 as Bounded>::min_value();
        if !max_cardinality {
            delta = self.compute_delta_vertices();
            self.update_mode = UpdateMode::Vertex;
        }
        let (delta, best_edge) =
            self.compute_delta_s_vertex_free_vertex(delta, &edges, &weight_matrix);
        let (delta, best_edge) =
            self.compute_delta_s_blossoms(delta, best_edge, edges, weight_matrix);
        let (update_blossom, mut delta) = self.compute_delta_t_blossoms(delta);

        if self.update_mode == UpdateMode::Undetermined {
            assert!(max_cardinality);
            // no more updates necessary
            self.update_mode = UpdateMode::Vertex;
            delta = if delta < 0.0 { 0.0 } else { delta };
        }
        dbg!(self.update_mode);

        (delta, best_edge, update_blossom)
    }

    pub fn update_blossom_structure(
        &mut self,
        stack: &mut Vec<usize>,
        best_edge: usize,
        update_blossom: usize,
        edges: &Vec<(usize, usize)>,
        endpoints: &Vec<usize>,
        matching: &Vec<i64>,
    ) -> bool {
        match self.update_mode {
            UpdateMode::Undetermined => (),
            UpdateMode::Vertex => return true,
            UpdateMode::SVertexFreeVertex => {
                self.allowed_edge[best_edge] = true;
                let (mut i, mut j) = edges[best_edge];
                if self.blossom_labels[self.blossom_id[i]] == 0 {
                    (j, i) = (i, j);
                }
                assert_eq!(self.blossom_labels[self.blossom_id[i]], 1);
                stack.push(i);
            }
            UpdateMode::SBlossom => {
                self.allowed_edge[best_edge] = true;
                let (i, j) = edges[best_edge];
                assert_eq!(self.blossom_labels[self.blossom_id[i]], 1);
                stack.push(i);
            }
            UpdateMode::TBlossom => {
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

    pub fn collect_leaves(&self, current_blossom: usize) -> Vec<usize> {
        // another heuristic to balance cache hits vs heap allocations
        let mut leaves = Vec::with_capacity(self.n_vertices >> 3);
        for v in BlossomLeaves::new(current_blossom, &self) {
            leaves.push(v);
        }
        leaves
    }

    pub fn assign_label<'a>(
        &'a mut self,
        stack: &mut Vec<usize>,
        vertex: usize,
        label: u8,
        endpoint: i64,
        matching: &Vec<i64>,
    ) {
        println!("assign_label: {} {} {}", vertex, label, endpoint);
        let b = self.blossom_id[vertex];
        assert_eq!(self.blossom_labels[vertex], 0);
        assert_eq!(self.blossom_labels[b], 0);
        self.blossom_labels[vertex] = label;
        self.blossom_labels[b] = label;
        self.label_endpoints[vertex] = endpoint;
        self.label_endpoints[b] = endpoint;
        self.best_edge[vertex] = -1;
        self.best_edge[b] = -1;
        if label == 1 {
            println!("Append: {:?}", self.collect_leaves(b));
            stack.append(&mut self.collect_leaves(b));
        } else if label == 2 {
            let root = self.blossom_root[b] as usize;
            assert!(matching[root] >= 0);
            self.assign_label(
                stack,
                self.label_endpoints[matching[root] as usize] as usize,
                1,
                matching[root] ^ 1,
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
        let mut v = v as i64;
        let mut w = w as i64;
        let mut path = Vec::new();
        while v != -1 || w != -1 {
            let mut b = self.blossom_id[v as usize];
            if self.blossom_labels[b] == 2 {
                for c in path.iter() {
                    self.blossom_labels[*c] = 1;
                }
                return Some(self.blossom_root[b] as usize);
            }
            assert_eq!(self.blossom_labels[b], 1);
            path.push(b);
            // why 5???
            self.blossom_labels[b] = 5;
            assert_eq!(
                self.label_endpoints[b],
                matching[self.blossom_root[b] as usize]
            );
            if self.label_endpoints[b] == -1 {
                v = -1;
            } else {
                v = endpoints[self.label_endpoints[b] as usize] as i64;
                b = self.blossom_id[v as usize];
                assert_eq!(self.blossom_labels[b], 2);
                assert!(self.label_endpoints[b] >= 0);
                v = endpoints[self.label_endpoints[b] as usize] as i64;
            }
            if w != -1 {
                let wtmp = w;
                w = v;
                v = wtmp;
            }
        }
        None
    }

    pub fn add_blossom<'a>(
        &'a mut self,
        stack: &mut Vec<usize>,
        root: usize,
        edge_idx: usize,
        edges: &Vec<(usize, usize)>,
        matching: &Vec<i64>,
        endpoints: &Vec<usize>,
        nbrhd_endpoints: &Vec<Vec<usize>>,
        weight_matrix: &SquareMatrix<f64>,
    ) {
        let (mut v, mut w) = edges[edge_idx];
        let bb = self.blossom_id[root];
        let mut bv = self.blossom_id[v];
        let mut bw = self.blossom_id[w];
        let b = match self.unused_blossoms.pop() {
            None => panic!("No more unused blossoms."),
            Some(c) => c,
        };

        self.blossom_root[b] = root as i64;
        self.blossom_parent[b] = -1;
        self.blossom_parent[bb] = b as i64;

        self.blossom_children[b] = Vec::new();
        self.blossom_endpoints[b] = Vec::new();

        while bv != bb {
            self.blossom_parent[bv] = b as i64;
            self.blossom_children[b].push(bv);
            self.blossom_endpoints[b].push(self.label_endpoints[bv]);
            assert!(
                self.blossom_labels[bv] == 2
                    || (self.blossom_labels[bv] == 1
                        && self.label_endpoints[bv] == matching[self.blossom_root[bv] as usize])
            );
            assert!(self.label_endpoints[bv] >= 0);
            v = endpoints[self.label_endpoints[bv] as usize];
            bv = self.blossom_id[v];
        }

        self.blossom_children[b].push(bb);
        self.blossom_children[b] = self.blossom_children[b].iter().copied().rev().collect();
        self.blossom_endpoints[b] = self.blossom_endpoints[b].iter().copied().rev().collect();
        self.blossom_endpoints[b].push(2 * edge_idx as i64);

        while bw != bb {
            self.blossom_parent[bw] = b as i64;
            self.blossom_children[b].push(bw);
            self.blossom_endpoints[b].push(self.label_endpoints[bw]);
            assert!(
                self.blossom_labels[bw] == 2
                    || (self.blossom_labels[bw] == 1
                        && self.label_endpoints[bw] == matching[self.blossom_root[bw] as usize])
            );
            assert!(self.label_endpoints[bw] >= 0);
            w = endpoints[self.label_endpoints[bw] as usize];
            bw = self.blossom_id[w];
        }

        assert_eq!(self.blossom_labels[bb], 1);
        self.blossom_labels[b] = 1;
        self.label_endpoints[b] = self.label_endpoints[bb];
        self.dual_soln[b] = 0.0;
        let leaves: Vec<usize> = self.collect_leaves(b);
        let mut t_leaves = Vec::with_capacity(leaves.len());
        for leaf in leaves.iter() {
            if self.blossom_id[*leaf] == 2 {
                t_leaves.push(*leaf);
            }
            self.blossom_id[*leaf] = b;
        }

        let mut best_edges = vec![-1i64; 2 * self.n_vertices];
        for bv in self.blossom_children[b].iter() {
            if self.blossom_best_edges[*bv].len() == 0 {
                for nbr_list in BlossomLeaves::new(*bv, &self)
                    .map(|v| nbrhd_endpoints[v].iter().map(|p| *p / 2))
                {
                    for nbr in nbr_list {
                        let (mut i, mut j) = edges[nbr];
                        if self.blossom_id[j] == b {
                            (j, i) = (i, j);
                        }
                        let bj = self.blossom_id[j];
                        if bj != b
                            && self.blossom_labels[bj] == 1
                            && (best_edges[bj] == -1
                                || self.slack(nbr, edges, weight_matrix)
                                    < self.slack(best_edges[bj] as usize, edges, weight_matrix))
                        {
                            best_edges[bj] = nbr as i64;
                        }
                    }
                }
            } else {
                for nbr in self.blossom_best_edges[*bv].iter() {
                    let (mut i, mut j) = edges[*nbr];
                    if self.blossom_id[j] == b {
                        (j, i) = (i, j);
                    }
                    let bj = self.blossom_id[j];
                    if bj != b
                        && self.blossom_labels[bj] == 1
                        && (best_edges[bj] == -1
                            || self.slack(*nbr, edges, weight_matrix)
                                < self.slack(best_edges[bj] as usize, edges, weight_matrix))
                    {
                        best_edges[bj] = *nbr as i64;
                    }
                }
            }
            self.blossom_best_edges[*bv] = Vec::new();
            self.best_edge[*bv] = -1;
        }
        self.blossom_best_edges[b] = best_edges
            .iter()
            .filter(|k| **k > -1)
            .map(|k| *k as usize)
            .collect();
        self.best_edge[b] = -1;
        for k in self.blossom_best_edges[b].iter() {
            if self.best_edge[b] == -1
                || self.slack(*k, edges, weight_matrix)
                    < self.slack(self.best_edge[b] as usize, edges, weight_matrix)
            {
                self.best_edge[b] = *k as i64;
            }
        }

        stack.append(&mut t_leaves);
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
                let leaves: Vec<usize> = self.collect_leaves(c);
                for v in leaves {
                    self.blossom_id[v] = c;
                }
            }
        }

        if !endstage && self.blossom_labels[blossom_id] == 2 {
            assert!(self.label_endpoints[blossom_id] >= 0);
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
                    self.assign_label(stack, endpoints[ix], 2, endpoint, matching);
                    self.allowed_edge[(self.blossom_endpoints[blossom_id][k] / 2) as usize] = true;
                    k += 1;
                    endpoint = self.blossom_endpoints[blossom_id][k];
                    self.allowed_edge[(endpoint / 2) as usize] = true;
                    k += 1;
                }

                k = 0;
                let mut bv = self.blossom_children[blossom_id][0];
                let ix = (endpoint ^ 1) as usize;
                self.blossom_labels[endpoints[ix]] = 2;
                self.blossom_labels[bv] = 2;
                self.label_endpoints[endpoints[ix]] = endpoint;
                self.label_endpoints[bv] = endpoint;
                self.best_edge[bv] = -1;

                // why not start at 0 here???
                // careful!
                k = 1;
                while self.blossom_children[blossom_id][j] != entry_child {
                    bv = self.blossom_children[blossom_id][j];
                    if self.blossom_labels[bv] == 1 {
                        k += 1;
                    } else {
                        let mut v = None;
                        for leaf in BlossomLeaves::new(bv, &self) {
                            if self.blossom_labels[leaf] != 0 {
                                v = Some(leaf);
                                break;
                            }
                        }
                        match v {
                            None => continue,
                            Some(v) => {
                                assert_eq!(self.blossom_labels[v], 2);
                                assert_eq!(self.blossom_id[v], bv);
                                self.blossom_labels[v] = 0;
                                self.blossom_labels[endpoints
                                    [matching[self.blossom_root[bv] as usize] as usize]] = 0;
                                self.assign_label(stack, v, 2, endpoint, matching);
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
                    self.assign_label(stack, endpoints[ix], 2, endpoint, matching);
                    self.allowed_edge[(self.blossom_endpoints[blossom_id][k - 1] / 2) as usize] =
                        true;
                    k -= 1;
                    endpoint = self.blossom_endpoints[blossom_id][k - 1];
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
                        for leaf in BlossomLeaves::new(bv, &self) {
                            if self.blossom_labels[leaf] != 0 {
                                v = Some(leaf);
                                break;
                            }
                        }
                        match v {
                            None => continue,
                            Some(v) => {
                                assert_eq!(self.blossom_labels[v], 2);
                                assert_eq!(self.blossom_id[v], bv);
                                self.blossom_labels[v] = 0;
                                self.blossom_labels[endpoints
                                    [matching[self.blossom_root[bv] as usize] as usize]] = 0;
                                self.assign_label(stack, v, 2, endpoint, matching);
                            }
                        }
                        k -= 1;
                    }
                }
            }
        }

        // recycle the blossom id
        self.blossom_labels[blossom_id] = 0;
        self.label_endpoints[blossom_id] = -1;
        self.blossom_root[blossom_id] = -1;
        self.blossom_best_edges[blossom_id] = Vec::new();
        self.best_edge[blossom_id] = -1;
        self.unused_blossoms.push(blossom_id);
    }

    // this function simultaneously expands blossoms
    // and appropriately augments the matching
    pub fn augment_blossom(
        &mut self,
        blossom_id: usize,
        vertex: usize,
        endpoints: &Vec<usize>,
        matching: &mut Vec<i64>,
    ) {
        // let t be the blossom which contains vertex
        // and which is a direct child of blossom_id
        let mut t = vertex;
        while self.blossom_parent[t] as usize != blossom_id {
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
        let mut j = i;
        if i % 2 == 1 {
            while j < self.blossom_children[blossom_id].len() {
                j += 1;
                t = self.blossom_children[blossom_id][j];
                let p = self.blossom_endpoints[blossom_id][j] as usize;
                let p_parity = p ^ 1;
                if t >= self.n_vertices {
                    self.augment_blossom(t, endpoints[p], endpoints, matching);
                }

                j += 1;
                t = self.blossom_children[blossom_id][j];
                if t >= self.n_vertices {
                    self.augment_blossom(t, endpoints[p_parity], endpoints, matching);
                }

                matching[endpoints[p]] = p_parity as i64;
                matching[endpoints[p_parity]] = p as i64;
            }
        } else {
            while j > 0 {
                j -= 1;
                t = self.blossom_children[blossom_id][j];
                let p = (self.blossom_endpoints[blossom_id][j - 1] ^ 1) as usize;
                let p_parity = p ^ 1;
                if t >= self.n_vertices {
                    self.augment_blossom(t, endpoints[p], endpoints, matching);
                }

                j -= 1;
                t = self.blossom_children[blossom_id][j];
                if t >= self.n_vertices {
                    self.augment_blossom(t, endpoints[p_parity], endpoints, matching);
                }

                matching[endpoints[p]] = p_parity as i64;
                matching[endpoints[p_parity]] = p as i64;
            }
        }

        // better to use BTReeSets?
        self.blossom_children[blossom_id].remove(i);
        self.blossom_endpoints[blossom_id].remove(i);
        self.blossom_root[blossom_id] = self.blossom_root[self.blossom_children[blossom_id][0]];
        assert_eq!(self.blossom_root[blossom_id] as usize, vertex);
    }

    pub fn augment_matching(
        &mut self,
        edge_idx: usize,
        edges: &Vec<(usize, usize)>,
        endpoints: &Vec<usize>,
        matching: &mut Vec<i64>,
    ) {
        let (v, w) = edges[edge_idx];
        let (mut s, mut p) = (v, 2 * edge_idx + 1);
        loop {
            let bs = self.blossom_id[s];
            assert_eq!(self.blossom_labels[bs], 1);
            assert_eq!(
                self.label_endpoints[bs],
                matching[self.blossom_root[bs] as usize]
            );
            if bs >= self.n_vertices {
                self.augment_blossom(bs, s, endpoints, matching);
            }
            matching[s] = p as i64;

            if self.label_endpoints[bs] == -1 {
                break;
            }

            let t = endpoints[self.label_endpoints[bs] as usize];
            let bt = self.blossom_id[t];
            assert_eq!(self.blossom_labels[bt], 2);
            assert!(self.label_endpoints[bt] >= 0);
            s = endpoints[self.label_endpoints[bt] as usize];
            let j = endpoints[(self.label_endpoints[bt] as usize) ^ 1];
            assert_eq!(self.blossom_root[bt] as usize, t);
            if bt >= self.n_vertices {
                self.augment_blossom(bt, j, endpoints, matching);
            }
            matching[j] = self.label_endpoints[bt];
            p = (self.label_endpoints[bt] ^ 1) as usize;
        }

        let (mut s, mut p) = (w, 2 * edge_idx);
        loop {
            let bs = self.blossom_id[s];
            assert_eq!(self.blossom_labels[bs], 1);
            assert_eq!(
                self.label_endpoints[bs],
                matching[self.blossom_root[bs] as usize]
            );
            if bs >= self.n_vertices {
                self.augment_blossom(bs, s, endpoints, matching);
            }
            matching[s] = p as i64;

            if self.label_endpoints[bs] == -1 {
                break;
            }

            let t = endpoints[self.label_endpoints[bs] as usize];
            let bt = self.blossom_id[t];
            assert_eq!(self.blossom_labels[bt], 2);
            assert!(self.label_endpoints[bt] >= 0);
            s = endpoints[self.label_endpoints[bt] as usize];
            let j = endpoints[(self.label_endpoints[bt] as usize) ^ 1];
            assert_eq!(self.blossom_root[bt] as usize, t);
            if bt >= self.n_vertices {
                self.augment_blossom(bt, j, endpoints, matching);
            }
            matching[j] = self.label_endpoints[bt];
            p = (self.label_endpoints[bt] ^ 1) as usize;
        }
    }
}

pub struct RootedSubGraphForest<N: Neighborhood + FromIterator<usize> + Clone> {
    pub roots: Vec<usize>,
    pub vertex_index_remapping: BTreeMap<usize, usize>,
    pub adjacencies: Vec<N>,
    pub root_map: Vec<usize>,
    pub root_paths: Vec<Vec<usize>>,
}

impl<N: Neighborhood + FromIterator<usize> + Clone> RootedSubGraphForest<N> {
    pub fn new(roots: Vec<usize>) -> Self {
        let mut vertex_index_remapping = BTreeMap::new();
        for (i, v) in roots.iter().enumerate() {
            vertex_index_remapping.insert(*v, i);
        }
        RootedSubGraphForest {
            roots: roots.clone(),
            vertex_index_remapping: vertex_index_remapping,
            adjacencies: vec![N::new(); roots.len()],
            root_map: roots.clone(),
            root_paths: roots.iter().map(|v| vec![*v]).collect(),
        }
    }

    pub fn get_neighborhood(&self, vertex: usize) -> Option<&N> {
        let index = self.vertex_index_remapping.get(&vertex)?;
        Some(&self.adjacencies[*index])
    }

    pub fn get_vertex_index(&self, vertex: usize) -> usize {
        match self.vertex_index_remapping.get(&vertex) {
            Some(i) => *i,
            None => panic!("Attempted to find non-existent vertex in subgraph."),
        }
    }

    pub fn get_root(&self, vertex: usize) -> Option<usize> {
        let index = self.vertex_index_remapping.get(&vertex)?;
        Some(self.root_map[*index])
    }

    pub fn get_root_path(&self, vertex: usize) -> Option<&Vec<usize>> {
        let index = self.vertex_index_remapping.get(&vertex)?;
        Some(&self.root_paths[*index])
    }

    pub fn contains_vertex(&self, vertex: usize) -> bool {
        self.vertex_index_remapping.contains_key(&vertex)
    }

    pub fn maybe_get_distance(&self, vertex: usize) -> Option<usize> {
        let index = self.vertex_index_remapping.get(&vertex)?;
        Some(self.root_paths[*index].len() - 1)
    }

    pub fn add_edge(&mut self, edge: (usize, usize)) {
        // assumes first vertex is already in the tree
        // and that the second vertex is not
        self.vertex_index_remapping
            .insert(edge.1, self.vertex_index_remapping.len());
        let index = match self.vertex_index_remapping.get(&edge.0) {
            None => panic!("Attempted to add an edge to RootedSubGraphForest whose first vertex was not already in the subgraph."),
            Some(i) => i
        };
        self.adjacencies[*index].insert(edge.1);
        let mut new_adjacency = N::new();
        new_adjacency.insert(edge.0);
        self.adjacencies.push(new_adjacency);
        self.root_map.push(self.root_map[*index]);
        let mut new_root_path = self.root_paths[*index].clone();
        new_root_path.push(edge.1);
        self.root_paths.push(new_root_path);
    }

    pub fn check_for_unmarked_vertex_w_even_dist(
        &self,
        marked_vertices: &HashSet<usize>,
    ) -> Option<usize> {
        // find a vertex in the forest which has an even distance from its root
        // and which is contained in `marked_vertices`
        for (vertex, i) in self.vertex_index_remapping.iter() {
            if self.root_paths[*i].len() % 2 == 1 && !marked_vertices.contains(vertex) {
                return Some(*vertex);
            }
        }
        None
    }
}
