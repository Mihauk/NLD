use criterion::{criterion_group, criterion_main, Criterion};

use nld::with_threads;

fn bench_with_threads(c: &mut Criterion) {
    c.bench_function("bench_with_threads", |b| {
        b.iter(|| {
            for _ in 1..=10 {
                std::hint::black_box(with_threads(10, 100, 999, 100, 0.1, 0.5));
            }
        });
    });
}

criterion_group!(
    benches,
    bench_with_threads,
);
criterion_main!(benches);