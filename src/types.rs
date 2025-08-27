use num::Bounded;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt::Display;

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

impl<N: Neighborhood + FromIterator<usize> + Clone> Graph<(), N> {
    pub fn from_edges(edges: Vec<(usize, usize)>, directed: bool) -> Self {
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
}

impl<
        W: Bounded + Sized + Copy,
        N: Neighborhood + FromIterator<usize> + Clone,
    > Graph<W, N>
{
    pub fn from_weighted_edges(edges: Vec<Edge<W>>, directed: bool) -> Self {
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

pub struct RootedSubGraphForest<
    N: Neighborhood + FromIterator<usize> + Clone,
> {
    pub roots: Vec<usize>,
    pub vertex_index_remapping: BTreeMap<usize, usize>,
    pub adjacencies: Vec<N>,
    pub root_map: Vec<usize>,
    pub root_paths: Vec<Vec<usize>>,
}

impl<N: Neighborhood + FromIterator<usize> + Clone>
    RootedSubGraphForest<N>
{
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
        self.vertex_index_remapping.insert(edge.1, self.vertex_index_remapping.len());
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
