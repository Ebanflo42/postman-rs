use postman::algorithms::*;

fn main() {
    let weighted_edges = vec![
        (1, 2, 2.),
        (1, 3, -2.),
        (2, 3, 1.),
        (2, 4, -1.),
        (3, 4, -6.),
    ];
    let matching = max_weight_matching(&weighted_edges, true);
}
