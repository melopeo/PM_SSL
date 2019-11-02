## Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs

MATLAB implementation of the paper:

[P. Mercado, F. Tudisco, and M. Hein, Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs. In NeurIPS 2019.](https://github.com/melopeo/PM_SSL/blob/master/PaperAndPoster/paper_Short.pdf)

## Content:
- `example.m` : contains an easy example showing how to use the code

- `run_everything.m` : runs experiments contained in our [paper](https://github.com/melopeo/PM_SSL/blob/master/PaperAndPoster/paper_Short.pdf)
 
## Usage:
Let `Wcell` be a cell with the adjacency matrices of each layer , `p` the power of the power mean Laplacian, `y` an array with the class of labeled nodes (zero denotes node is unlabeled). Classes through the power mean Laplacian `L_p` regularizer are computed via
```
y_hat           = SSL_multilayer_graphs_with_power_mean_laplacian(Wcell, p, y);
```

## Quick Overview:
![](https://github.com/melopeo/PM_SSL/blob/master/PaperAndPoster/poster.jpg)

## Citation:
```
@article{mercadoNeurIPS2019,
  title = 	 {Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs},
  author = 	 {Mercado, Pedro and Tudisco, Francesco and Hein, Matthias},
  conference = 	 {NeurIPS},
  year = 	 {2019},
}

```
