use core::f64;
use std::collections::BTreeSet;

pub type EdgeList = Vec<(usize, usize)>;
pub type AdjacencyList = Vec<Vec<usize>>;
pub type WeightedEdgeList = Vec<(usize, usize, f64)>;
pub type WeightedAdjacencyList = Vec<Vec<(usize, f64)>>;

pub trait Container {
    type Item;
    fn contains(&self, item: Self::Item) -> bool;
}

pub struct Graph {
    vertices: Vec<usize>,
    adjacency_list: Vec<BTreeSet<usize>>,
}

pub struct WeightedGraph {
    vertices: Vec<usize>,
    adjacency_list: Vec<Vec<(usize, f64)>>,
}

pub struct WeightMatrix {
    data: Vec<f64>,
    num_vertices: usize,
}

fn get_num_vertices(edge_list: &EdgeList) -> usize {
    let mut result = 0;
    for edge in edge_list.iter() {
        if edge.0 > result {
            result = edge.0;
        }
        if edge.1 > result {
            result = edge.1;
        }
    }
    result + 1
}

fn get_num_vertices_weighted(weighted_edge_list: &WeightedEdgeList) -> usize {
    let mut result = 0;
    for edge in weighted_edge_list.iter() {
        if edge.0 > result {
            result = edge.0;
        }
        if edge.1 > result {
            result = edge.1;
        }
    }
    result + 1
}

impl WeightMatrix {
    pub fn new(weighted_edge_list: WeightedEdgeList, directed: bool, no_edge_value: f64) -> Self {
        if weighted_edge_list.len() == 0 {
            return WeightMatrix {
                data: Vec::new(),
                num_vertices: 0,
            };
        }

        let num_vertices = get_num_vertices_weighted(&weighted_edge_list);

        let mut result_data = vec![no_edge_value; num_vertices * num_vertices];
        for (i, j, w) in weighted_edge_list.iter() {
            result_data[num_vertices * i + j] = *w;
            if !directed {
                result_data[num_vertices * j + i] = *w;
            }
        }

        WeightMatrix {
            data: result_data,
            num_vertices: num_vertices,
        }
    }

    pub fn get_ix(&self, i: usize, j: usize) -> f64 {
        self.data[self.num_vertices * i + j]
    }

    pub fn set_ix(&mut self, i: usize, j: usize, val: f64) {
        self.data[self.num_vertices * i + j] = val;
    }

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
            .chunks(self.num_vertices)
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

pub fn adjacencies_from_edges(edge_list: EdgeList, directed: bool) -> AdjacencyList {
    let num_vertices = get_num_vertices(&edge_list);
    let mut result = Vec::with_capacity(num_vertices);
    for _ in 0..num_vertices {
        result.push(Vec::new());
    }

    for (i, j) in edge_list.iter() {
        result[*i].push(*j);
        if !directed {
            result[*j].push(*i);
        }
    }

    result
}

pub fn weighted_adjacencies_from_edges(
    weighted_edge_list: WeightedEdgeList,
    directed: bool,
) -> WeightedAdjacencyList {
    let num_vertices = get_num_vertices_weighted(&weighted_edge_list);
    let mut result = Vec::with_capacity(num_vertices);
    for _ in 0..num_vertices {
        result.push(Vec::new());
    }

    for (i, j, w) in weighted_edge_list.iter() {
        result[*i].push((*j, *w));
        if !directed {
            result[*j].push((*i, *w));
        }
    }

    result
}

pub fn floyd_warshall(weighted_edge_list: WeightedEdgeList, directed: bool) -> WeightMatrix {
    let mut result = WeightMatrix::new(weighted_edge_list, directed, f64::INFINITY);
    println!("{}", result.to_string());
    for i in 0..result.num_vertices {
        result.set_ix(i, i, 0f64);
    }

    //let mut pathtracker = Vec::with_capacity(result.num_vertices*result.num_vertices);
    //for i in 0..result.num_vertices {
    //    pathtracker.extend(vec![i; result.num_vertices])
    //}

    for j in 0..result.num_vertices {
        for i in 0..result.num_vertices {
            if i != j {
                for k in 0..result.num_vertices {
                    let sum = result.get_ix(i, j) + result.get_ix(j, k);
                    if k != j && result.get_ix(i, k) > sum {
                        result.set_ix(i, k, sum);
                        //pathtracker[result.num_vertices*i + k] = pathtracker[result.num_vertices*j + k];
                    }
                }
            }
        }
    }

    result
}

struct RootedSubGraphForest {
    roots: Vec<usize>,
    vertex_index_remapping: Vec<usize>,
    adjacencies: AdjacencyList,
    root_map: Vec<usize>,
    root_paths: Vec<Vec<usize>>,
}

impl RootedSubGraphForest {
    fn new(roots: &Vec<usize>) -> Self {
        RootedSubGraphForest {
            roots: roots.clone(),
            vertex_index_remapping: roots.clone(),
            adjacencies: vec![Vec::new(); roots.len()],
            root_map: roots.clone(),
            root_paths: roots.iter().map(|v| vec![*v]).collect(),
        }
    }

    fn get_neighborhood(&self, vertex: usize) -> Option<&Vec<usize>> {
        let adjacency_ix = self
            .vertex_index_remapping
            .iter()
            .position(|v| *v == vertex)?;
        Some(&self.adjacencies[adjacency_ix])
    }

    fn get_vertex_index(&self, vertex: usize) -> usize {
        match self
            .vertex_index_remapping
            .iter()
            .position(|v| *v == vertex)
        {
            Some(i) => i,
            None => panic!("Attempted to find non-existent vertex in subgraph."),
        }
    }

    fn get_root(&self, vertex: usize) -> Option<usize> {
        let adjacency_ix = self
            .vertex_index_remapping
            .iter()
            .position(|v| *v == vertex)?;
        Some(self.root_map[adjacency_ix])
    }

    fn get_root_path(&self, vertex: usize) -> Option<Vec<usize>> {
        let adjacency_ix = self
            .vertex_index_remapping
            .iter()
            .position(|v| *v == vertex)?;
        Some(self.root_paths[adjacency_ix].clone())
    }

    fn contains_vertex(&self, vertex: usize) -> bool {
        self.vertex_index_remapping.iter().any(|v| *v == vertex)
    }

    fn maybe_get_distance(&self, vertex: usize) -> Option<usize> {
        let adjacency_ix = self
            .vertex_index_remapping
            .iter()
            .position(|v| *v == vertex)?;
        Some(self.root_paths[adjacency_ix].len() - 1)
    }

    fn add_edge(&mut self, edge: (usize, usize)) {
        // assumes first vertex is already in the tree
        // and that the second index is not
        self.vertex_index_remapping.push(edge.1);
        let adjacency_ix = match self.vertex_index_remapping.iter().position(|v| *v == edge.0) {
            Some(i) => i,
            None => panic!("Attempted to add an edge to a rooted tree where the first vertex of the edge was not already in the tree")
        };
        self.adjacencies[adjacency_ix].push(edge.1);
        self.adjacencies.push(vec![edge.0]);
        self.root_map.push(self.root_map[edge.0]);
        let mut new_root_path = self.root_paths[edge.0].clone();
        new_root_path.push(edge.1);
        self.root_paths.push(new_root_path);
    }

    fn check_for_unmarked_edge(&self, marked_edges: &EdgeList, incident: usize) -> Option<usize> {
        // return an unmarked edge, if one exists
        let adjacency_ix = match self
            .vertex_index_remapping
            .iter()
            .position(|v| *v == incident)
        {
            Some(i) => i,
            None => {
                panic!("Attempted to query neighborhood of a vertex which is not in the subgraph")
            }
        };
        for v in self.adjacencies[adjacency_ix].iter() {
            if marked_edges
                .iter()
                .all(|e| *e != (*v, incident) && *e != (incident, *v))
            {
                return Some(*v);
            }
        }
        None
    }

    fn check_for_unmarked_vertex_w_even_dist(&self, marked_vertices: &Vec<usize>) -> Option<usize> {
        for (i, vertex) in self.vertex_index_remapping.iter().enumerate() {
            if self.root_paths[i].len() % 2 == 1 && marked_vertices.iter().all(|v| *v != *vertex) {
                return Some(*vertex);
            }
        }
        None
    }
}

fn contract_blossom(
    graph: &Graph,
    matching: &EdgeList,
    blossom: Vec<usize>,
) -> (Graph, EdgeList, usize) {
    // assume the first vertex in the blossom is unmatched
    // the blossom passed to this function should always satisfy this assumption
    let blossom_root = blossom[0];

    let mut new_vertices = graph.vertices.clone();
    for v in blossom.iter() {
        if *v != blossom_root {
            if let Some(i) = new_vertices.iter().position(|w| *w == *v) {
                new_vertices.remove(i);
            }
        }
    }

    let mut new_adjacencies = graph.adjacency_list.clone();
    for i in 0..new_adjacencies.len() {
        if i == blossom_root {
            for b in blossom.iter() {
                for v in graph.adjacency_list[*b].iter() {
                    new_adjacencies[i].insert(*v);
                }
            }
        } else if blossom.iter().any(|b| *b == i) {
            new_adjacencies[i].clear();
        } else {
            for b in blossom.iter() {
                if *b != blossom_root {
                    if new_adjacencies[i].remove(b) {
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
        },
        new_matching,
        blossom_root,
    )
}

fn expand_blossom(
    original_graph: &Graph,
    original_matching: &EdgeList,
    contracted_matching: &EdgeList,
    path: Vec<usize>,
    blossom: &Vec<usize>,
    blossom_root: usize,
) -> Vec<usize> {
    let blossom_contains = |v: usize| blossom.iter().any(|u| *u == v);

    if let Some(i) = path.iter().position(|v| *v == blossom_root) {
        if i == 0 {
            path
        } else if i == path.len() - 1 {
            path
        } else {
            // if the path passes through the blossom, there are two ways to expand the path
            // and we have to carefully choose the one which ensures that the resulting path
            // is still alternating
            let mut first_segment = Vec::from(&path[..i]);
            let incoming_vert = *first_segment.last().unwrap();
            let second_segment = Vec::from(&path[(i + 1)..]);
            let outgoing_vert = second_segment[0];

            let incoming_edge_matched = contracted_matching.iter().any(|e| {
                *e == (incoming_vert, blossom_root) || *e == (blossom_root, incoming_vert)
            });

            if incoming_edge_matched {
                let outgoing_blossom_vert_idx =
                    match blossom.into_iter().position(|v| original_graph.adjacency_list[outgoing_vert].contains(v)) {
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
                            == blossom
                                [(outgoing_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
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
                            == blossom
                                [(outgoing_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
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
                    match blossom.into_iter().position(|v| original_graph.adjacency_list[incoming_vert].contains(v)) {
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
                            == blossom
                                [(incoming_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
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
                            == blossom
                                [(incoming_blossom_vert_idx + blossom.len() - 1) % blossom.len()]
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
        }
    } else {
        path
    }
}

fn find_augmenting_path(
    graph: &Graph,
    matching: &EdgeList,
    exposed_vertices: &Vec<usize>,
) -> Vec<usize> {
    let mut marked_vertices = Vec::new();
    let mut marked_edges = matching.clone();
    let mut forest = RootedSubGraphForest::new(exposed_vertices);

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
                            let (contracted_graph, contracted_matching, blossom_idx) =
                                contract_blossom(graph, matching, vw_path);
                            let contracted_path = find_augmenting_path(
                                &contracted_graph,
                                &contracted_matching,
                                exposed_vertices,
                            );
                        }
                    }
                }
            }
            marked_edges.push((v, w));
        }
        marked_vertices.push(v);
    }
    Vec::new()
}

pub fn edmonds_max_cardinality_matching(graph: &Graph) -> EdgeList {
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
                exposed_vertices.remove(i);
            } else {
                panic!("First vertex of augmenting path was not an exposed vertex!");
            }
            if let Some(j) = exposed_vertices.iter().position(|u| *u == w) {
                exposed_vertices.remove(j);
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
