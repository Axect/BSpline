use peroxide::fuga::*;
use anyhow::Result;

fn main() -> Result<()> {
    let knots = vec![0f64, 0f64, 0f64, 0f64, 1f64, 2f64, 3f64, 3f64, 3f64, 3f64];
    let degree = 3;
    let num_points = 100;
    let control_points = vec![
        vec![0f64, 2f64],
        vec![0.2, -1f64],
        vec![0.4, 1f64],
        vec![0.6, -1f64],
        vec![0.8, 1f64],
        vec![1f64, 2f64],
    ];

    let spline = BSpline::new(degree, knots, num_points, control_points);
    let bases = spline.generate_bases();
    spline.bases_plot(bases)?;
    spline.plot()?;

    Ok(())
}

struct BSpline {
    degree: usize,
    knots: Vec<f64>,
    t_values: Vec<f64>,
    control_points: Vec<Vec<f64>>,
}

impl BSpline {
    fn new(degree: usize, knots: Vec<f64>, num_points: usize, control_points: Vec<Vec<f64>>) -> Self {
        let t_values = linspace(knots[0], knots[knots.len()-1], num_points);
        assert_eq!(control_points.len(), knots.len() - (degree + 1));
        BSpline { degree, knots, t_values, control_points }
    }

    //fn cox_de_boor(&self, t: f64, p: usize, i: usize) -> f64 {
    //    if p == 0 {
    //        if (self.knots[i] <= t && t < self.knots[i + 1]) || (i == self.knots.len() - (self.degree + 2) && t == self.knots[i + 1]) {
    //            1f64
    //        } else {
    //            0f64
    //        }
    //    } else {
    //        let a = if self.knots[i + p] == self.knots[i] {
    //            0f64
    //        } else {
    //            (t - self.knots[i]) / (self.knots[i + p] - self.knots[i])
    //        };
    //        let b = if self.knots[i + p + 1] == self.knots[i + 1] {
    //            0f64
    //        } else {
    //            (self.knots[i + p + 1] - t) / (self.knots[i + p + 1] - self.knots[i + 1])
    //        };
    //        a * self.cox_de_boor(t, p - 1, i) + b * self.cox_de_boor(t, p - 1, i + 1)
    //    }
    //}
    fn cox_de_boor(&self, t: f64, p: usize, i: usize) -> f64 {
        let mut n = vec![vec![0.0; p + 1]; p + 1];

        // Initialize the zeroth degree basis functions
        for (j, n_j) in n.iter_mut().enumerate() {
            if (self.knots[i + j] <= t && t < self.knots[i + j + 1]) || (i + j == self.knots.len() - (self.degree + 2) && t == self.knots[i + j + 1]) {
                n_j[0] = 1.0;
            } else {
                n_j[0] = 0.0;
            }
        }

        // Compute the basis functions for higher degrees
        for k in 1..=p {
            for j in 0..=(p - k) {
                let a = if self.knots[i + j + k] == self.knots[i + j] {
                    0.0
                } else {
                    (t - self.knots[i + j]) / (self.knots[i + j + k] - self.knots[i + j])
                };

                let b = if self.knots[i + j + k + 1] == self.knots[i + j + 1] {
                    0.0
                } else {
                    (self.knots[i + j + k + 1] - t) / (self.knots[i + j + k + 1] - self.knots[i + j + 1])
                };

                n[j][k] = a * n[j][k - 1] + b * n[j + 1][k - 1];
            }
        }

        n[0][p]
    }

    fn generate_bases(&self) -> Vec<Vec<f64>> {
        let mut bases = vec![];
        for i in 0..self.knots.len() - (self.degree + 1) {
            bases.push(self.t_values.fmap(|x| self.cox_de_boor(x, self.degree, i)));
        }
        bases
    }

    fn bases_plot(&self, bases: Vec<Vec<f64>>) -> Result<()> {
        let mut plt = Plot2D::new();
        plt.set_domain(self.t_values.clone());
        for basis in bases {
            plt.insert_image(basis);
        }
        plt.set_xlabel("t")
            .set_ylabel("Basis")
            .set_style(PlotStyle::Nature)
            .tight_layout()
            .set_dpi(600)
            .set_path("plot.png")
            .savefig()?;
        Ok(())
    }

    fn interpolate(&self) -> (Vec<f64>, Vec<f64>) {
        let mut x = vec![];
        let mut y = vec![];
        for t in self.t_values.clone() {
            let mut x_val = 0f64;
            let mut y_val = 0f64;
            for i in 0..self.control_points.len() {
                let basis = self.cox_de_boor(t, self.degree, i);
                x_val += basis * self.control_points[i][0];
                y_val += basis * self.control_points[i][1];
            }
            x.push(x_val);
            y.push(y_val);
        }
        (x, y)
    }

    fn plot(&self) -> Result<()> {
        let (x, y) = self.interpolate();
        let mut plt = Plot2D::new();
        plt.set_domain(x)
            .insert_image(y)
            .set_xlabel(r"$x$")
            .set_ylabel(r"$y$")
            .set_style(PlotStyle::Nature)
            .tight_layout()
            .set_dpi(600)
            .set_path("interpolate_plot.png")
            .savefig()?;

        Ok(())
    }
}

