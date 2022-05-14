# QuasiBootsMedian
Quasi-boots meidan method for MR analyses

Please source the script median_boots.R in R, and use function: median_boots(bx,bxse,by,byse,cormat,same_population, IV_list, boots_iteration = 10000, seed=314159265) to run quasi-boots median method.

Input:
bx: genetic association with exposure
bxse: standard error of genetic association with exposure
by:genetic association with outcome
byse:standard error of genetic association with outcome
cormat:genetic correlation matrix, please add colnames and rownames as genetic variants names for this matrix (matching IV_list)
IV_list: genetic variants names (make sure the order for bx, bxse, by and byse matches this list)
same_population: True if one-sample MR; False if two-sample MR
boots_iteration: number of iterations for the bootstrap, default is 10000
seed: random seed, default is 314159265
