use core::f64;

pub type EdgeList = Vec<(usize, usize)>;
pub type AdjacencyList = Vec<Vec<usize>>;
pub type WeightedEdgeList = Vec<(usize, usize, f64)>;
pub type WeightedAdjacencyList = Vec<Vec<(usize, f64)>>;

pub struct WeightMatrix {
    data: Vec<f64>,
    num_vertices: usize
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
                num_vertices: 0
            };
        }

        let num_vertices = get_num_vertices_weighted(&weighted_edge_list);

        let mut result_data = vec![no_edge_value; num_vertices*num_vertices];
        for (i, j, w) in weighted_edge_list.iter() {
            result_data[num_vertices*i + j] = *w;
            if !directed {
                result_data[num_vertices*j + i] = *w;
            }
        }

        WeightMatrix {
            data: result_data,
            num_vertices: num_vertices
        }
    }

    pub fn get_ix(&self, i: usize, j: usize) -> f64 {
        self.data[self.num_vertices*i + j]
    }

    pub fn set_ix(&mut self, i: usize, j: usize, val: f64) {
        self.data[self.num_vertices*i + j] = val;
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
        let helper: Vec<String> = self.data
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

pub fn weighted_adjacencies_from_edges(weighted_edge_list: WeightedEdgeList, directed: bool) -> WeightedAdjacencyList {
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
            (0, 3, 1.0)
        ];
        let distance_matrix = floyd_warshall(weighted_edges, false);
        println!("{}", distance_matrix.to_string());
        assert_eq!(distance_matrix.data, vec![0.0, 1.0, 2.0, 1.0, 1.0,
                                              1.0, 0.0, 1.0, 2.0, 2.0,
                                              2.0, 1.0, 0.0, 2.0, 3.0,
                                              1.0, 2.0, 2.0, 0.0, 1.0,
                                              1.0, 2.0, 3.0, 1.0, 0.0]);
    }
}
