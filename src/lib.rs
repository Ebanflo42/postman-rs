pub mod algorithms;
pub mod types;

#[cfg(test)]
mod tests {
    use crate::algorithms::*;
    use crate::types::*;
    use std::collections::BTreeSet;
    use std::collections::HashSet;

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
    /*
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
            (16, 17),
        ];
        let graph: Graph<(), Vec<usize>> = Graph::from_edges(&edges, false);
        let matching = max_cardinality_matching(&graph);
        let true_matching: Vec<(usize, usize)> = vec![
            (0, 1),
            (2, 3),
            (8, 5),
            (6, 9),
            (4, 7),
            (10, 11),
            (12, 13),
            (14, 15),
            (16, 17),
        ];
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
        let graph: Graph<(), BTreeSet<usize>> =
            Graph::from_edges_and_vertices(&edges, &vertices, false);
        let matching = max_cardinality_matching(&graph);
        let true_matching: Vec<(usize, usize)> = vec![
            (0, 1),
            (2, 3),
            (8, 5),
            (6, 9),
            (4, 7),
            (10, 11),
            (12, 13),
            (14, 15),
        ];
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
            (16, 17),
        ];
        let vertices = (0usize..18).collect();
        let graph: Graph<(), HashSet<usize>> =
            Graph::from_edges_and_vertices(&edges, &vertices, false);
        let matching = max_cardinality_matching(&graph);
        let true_matching: Vec<(usize, usize)> = vec![
            (0, 1),
            (2, 3),
            (8, 5),
            (6, 9),
            (4, 7),
            (10, 11),
            (12, 13),
            (14, 15),
            (16, 17),
        ];
        println!("{:?}", matching);
        assert_eq!(matching, true_matching);
    }
    */
    #[test]
    fn max_weight_matching1() {
        let weighted_edges = vec![(1, 2, 5.), (2, 3, 11.), (3, 4, 5.)];
        let matching = max_weight_matching(&weighted_edges, false);
        assert_eq!(matching, vec![-1, -1, 3, 2, -1]);
        println!("{:?}", matching);
    }

    #[test]
    fn max_weight_matching2() {
        let weighted_edges = vec![
            (1, 2, 9.),
            (1, 3, 9.),
            (2, 3, 10.),
            (2, 4, 8.),
            (3, 5, 8.),
            (4, 5, 10.),
            (5, 6, 6.),
        ];
        let matching = max_weight_matching(&weighted_edges, false);
        assert_eq!(matching, vec![-1, 3, 4, 1, 2, 6, 5]);
        println!("{:?}", matching);
    }

    #[test]
    fn max_weight_matching3() {
        let weighted_edges = vec![
            (1, 2, 40.),
            (1, 3, 40.),
            (2, 3, 60.),
            (2, 4, 55.),
            (3, 5, 55.),
            (4, 5, 50.),
            (1, 8, 15.),
            (5, 7, 30.),
            (7, 6, 10.),
            (8, 10, 10.),
            (4, 9, 30.),
        ];
        let matching = max_weight_matching(&weighted_edges, false);
        assert_eq!(matching, vec![-1, 2, 1, 5, 9, 3, 7, 6, 10, 4, 8]);
        println!("{:?}", matching);
    }

    #[test]
    fn max_weight_matching4() {
        let weighted_edges = vec![
            (1, 2, 45.),
            (1, 7, 45.),
            (2, 3, 50.),
            (3, 4, 45.),
            (4, 5, 95.),
            (4, 6, 94.),
            (5, 6, 94.),
            (6, 7, 50.),
            (1, 8, 30.),
            (3, 11, 35.),
            (5, 9, 36.),
            (7, 10, 26.),
            (11, 12, 5.),
        ];
        let matching = max_weight_matching(&weighted_edges, false);
        //assert_eq!(matching, vec![-1, 2, 1, 5, 9, 3, 7, 6, 10, 4, 8]);
        println!("{:?}", matching);
    }
}
