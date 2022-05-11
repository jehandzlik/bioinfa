## 10-May-2022
## Normalization

## Data

## A) Same library size, different genes

# control genes
control = c(10,300,500,20,100,8)

# treatment genes
# two genes are de, one expressed more and one expressed less
treatment = c(10,300*20,500,20,100*0.1,8)

## B) Different library size, different genes
treatment = c(10,300*20,500,20,100*0.1,8) * 20


