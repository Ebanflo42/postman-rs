pub mod types;
pub mod algorithms;

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;
    use std::collections::HashSet;
    use crate::types::*;
    use crate::algorithms::*;

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
        let distance_matrix = floyd_warshall(&weighted_edges, false);
        println!("{}", distance_matrix.to_string());
        assert_eq!(
            distance_matrix.data,
            vec![
                0.0, 1.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 2.0, 3.0, 1.0,
                2.0, 2.0, 0.0, 1.0, 1.0, 2.0, 3.0, 1.0, 0.0
            ]
        );
    }

    #[test]
    fn edmonds1() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (4, 6),
            (4, 7),
            (5, 8),
            (6, 9),
            (7, 10),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 15),
            (15, 11),
            (15, 16),
            (16, 17)
        ];
        let graph: Graph<(), Vec<usize>> = Graph::from_edges(&edges, false);
        let matching = max_cardinality_matching(&graph);
        let true_matching: Vec<(usize, usize)> = vec![(0, 1), (2, 3), (8, 5), (6, 9), (4, 7), (10, 11), (12, 13), (14, 15), (16, 17)];
        println!("{:?}", matching);
        assert_eq!(matching, true_matching);
    }

    #[test]
    fn edmonds2() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (4, 6),
            (4, 7),
            (5, 8),
            (6, 9),
            (7, 10),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 15),
            (15, 11),
        ];
        let vertices = (0usize..17).collect();
        let graph: Graph<(), BTreeSet<usize>> = Graph::from_edges_and_vertices(&edges, &vertices, false);
        let matching = max_cardinality_matching(&graph);
        let true_matching: Vec<(usize, usize)> = vec![(0, 1), (2, 3), (8, 5), (6, 9), (4, 7), (10, 11), (12, 13), (14, 15)];
        println!("{:?}", matching);
        assert_eq!(matching, true_matching);
    }

    #[test]
    fn edmonds3() {
        let edges = vec![
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 15),
            (15, 11),
            (15, 16),
            (16, 17)
        ];
        let vertices = (0usize..18).collect();
        let graph: Graph<(), HashSet<usize>> = Graph::from_edges_and_vertices(&edges, &vertices, false);
        let matching = max_cardinality_matching(&graph);
        //let true_matching: Vec<(usize, usize)> = vec![(0, 1), (2, 3), (8, 5), (6, 9), (4, 7), (10, 11), (12, 13), (14, 15), (16, 17)];
        println!("{:?}", matching);
        //assert_eq!(matching, true_matching);
    }

}
