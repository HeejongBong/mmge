    //
    //  main.cpp
    //  matrix_grouping
    //
    //  Created by Yunzhang Zhu on 7/19/12.
    //  Copyright (c) 2012 Yunzhang Zhu. All rights reserved.
    //

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <R.h>
#include <Rinternals.h>

using namespace std;
#include "dmatrix.h"
#include "def.h"

extern "C"{
        // for debugging //
    void print_dmatrix(double* matrix,int m,int n){
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                printf("%5.3f ", matrix[j*m+i]);
            }
            printf("\n");
        }
        printf("\n");
    }
    
    void print_imatrix(int* matrix,int m,int n){
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                printf("%d ", matrix[j*m+i]);
            }
            printf("\n");
        }
        printf("\n");
    }
    
        // screening method with different lambda1's for each group//
    void screen_complete_new(int *AA, double* S_bar,double* lambda1, double* lambda2, double* nn, int* pp, int* LL){
        double* S = S_bar;
        int* A = AA;
        double lam1 = lambda1[0];
        double lam2 = lambda2[0];
        int p = pp[0];
        int L = LL[0];
        double* n = nn;
        double tmp1, tmp2;
        vector<vector<int> > subset;
        vector<int> tmpset;
        int tmp[L];
        for (int i = 1; i < 16; i++) {
            int j = i;
            for (int l = 0; l < L; l++) {
                tmp[l] = j%2;
                j = j/2;
            }
            tmpset.clear();
            for (int l = 0; l < L; l++) {
                if (tmp[l] == 1) tmpset.push_back(l);
            }
            subset.push_back(tmpset);
        }
        
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                double s[L];
                for (int l = 0; l < L; l++) {
                    s[l] = S[l*p*p+j*p+i]*n[l];
                }
                for (int k = 0; k < subset.size(); k++) {
                    tmpset = subset[k];
                    tmp1 = .0;
                    tmp2 = .0;
                    for (int l = 0; l < tmpset.size(); l++) {
                        tmp1 += s[tmpset[l]];
                        tmp2 += n[tmpset[l]];
                    }
                    if (abs(tmp1) > (lam1 * tmp2 + tmpset.size()*(L-tmpset.size())*lam2)) {
                        A[i*p+j] = 1;
                        break;
                    }
                }
            }
        }
    }
    
        // screening methods for complete graph situation //
    void screen_complete(int *AA, double* S_bar, double* lambda1, double* lambda2, double* nn, int* pp, int* LL){
        
        double* S = S_bar;
        int* A = AA;
        double lam1 = lambda1[0];
        double lam2 = lambda2[0];
        int p = pp[0];
        int L = LL[0];
        double* n = nn;
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                double s[L];
                for (int l = 0; l < L; l++) {
                    s[l] = S[l*p*p+j*p+i]*n[l];
                }
                sort(s, s + L);
                double sum1 = 0;
                double sum2  = 0;
                for (int l = 0; l < L; l++) {
                    sum1 += s[l];
                    sum2 += s[L-l-1];
                    double tmp1 = sum1/(l+1);
                    double tmp2 = sum2/(l+1);
                    if (max(abs(tmp1),abs(tmp2)) > (lam1 + (L-l-1)*lam2)){
                        A[i*p+j] = 1;
                        break;
                    }
                }
            }
        }
    }
    
    
        // solve .5*mu*|| x ||^2 - c^T x  + \sum_i lambda_i |x_i|
    void soft_thred(int n, double* x, double* c, double* mu, double* lambda){
        for (int i = 0; i < n; i++) {
            x[i] = sign(c[i])* max((abs(c[i]) - lambda[i]), 0.0) / mu[0];
        }
    }
    
        // check two vectors x and y if their non-zero patterns are the same;
        // return 1 if zero patterns are the same. //
    int zero_agree(int n, double* x, double* y){
        for (int i = 0; i < n; i++) {
            if ((x[i]==0) && (y[i]!=0)) {
                return 0;
            }
            
            if ((x[i]!=0) && (y[i]==0)) {
                return 0;
            }
        }
        return 1;
    }
    
        // updating matrix A based on covmat
	void update_A(double covmat[], double A[], int col, int p, int L){
            // cout << "updating matrix A based on covmat and covmat_inverse" << endl;
		for (int l = 0; l < L; l++) {
			double xjjl_inverse = 1/covmat[l*p*p+col*p+col];
			for (int i = 0; i < p-1; i++) {
				int i_index = ((i >= col)? (i+1):i);
				for (int j = 0; j < p-1; j++) {
					int j_index = ((j >= col)? (j+1):j);
					A[l*(p-1)*(p-1) + i*(p-1) + j] = covmat[l*p*p+i_index*p+j_index] - xjjl_inverse*covmat[l*p*p+col*p+j_index]*covmat[l*p*p+col*p+i_index];
				}
			}
		}
	}
    
        //updating covmat matrix based on covmat_inverse and A
	void update_cov(double covmat[], double covmat_inverse[], double A[], int col, int p, int L){
            //updating covmat matrix based on covmat_inverse and A
		for (int l = 0; l<L; l++) {
			double xjjl = covmat_inverse[l*p*p+col*p+col];
			
			double temp_a[p-1];
			for (int i = 0 ; i<p-1; i++) {
				temp_a[i] = 0;
				for (int j = 0; j<p-1; j++) {
					int j_index = ((j >= col)? (j+1):j);
					temp_a[i] += A[l*(p-1)*(p-1)+j*(p-1)+i]*covmat_inverse[l*p*p+col*p+j_index];
				}
			}
            
			for (int i = 0; i<p-1; i++) {
				int i_index = ((i >= col)? (i+1):i);
				xjjl -= covmat_inverse[l*p*p+col*p+i_index]*temp_a[i];
			}
			double xjjl_inverse = 1/xjjl;
			covmat[l*p*p+col*p+col] = xjjl_inverse;
			
			for (int i = 0; i<p-1; i++) {
				int i_index = ((i >= col)? (i+1):i);
				double temp = -temp_a[i]*xjjl_inverse;
				covmat[l*p*p+col*p+i_index] = temp;
				covmat[l*p*p+i_index*p+col] = temp;
				for (int j = 0; j<p-1; j++) {
					int j_index = ((j >= col)? (j+1):j);
					covmat[l*p*p+i_index*p+j_index] = A[l*(p-1)*(p-1)+i*(p-1)+j]+temp_a[j]*temp_a[i]*xjjl_inverse;
				}
			}
		}
	}
    
        //function evaluating \sum_l (1/2 a_l x_l^T A_l x_l + b_l^T x_l + lambda1 * sum_j |x_lj|) + lambda2 * sum_{l l' in graph} |x_l - x_l'|
	double group_fun(double *X, double *a,double *b,double *A,int *graph,double *lambda1_vec,double *lambda2_vec,int p,int L,int NumEdges){
		
		double part1 = 0; // l1 peanlty
		double part2 = 0; // grouping penalty
		double part3 = 0; // loss function
        double tmp[p-1];
        
        for (int l = 0; l < L; l++) {
            dmat_elemprod(p-1, X+l*(p-1), lambda1_vec+l*(p-1), tmp);
            part1 += dmat_norm1(p-1, tmp);
        }
        
		for (int m = 0; m < NumEdges; m++) {
            dmat_waxpby(p-1, 1.0, X+graph[2*m]*(p-1), -1.0, X+graph[2*m+1]*(p-1), tmp);
            dmat_elemprod(p-1, lambda2_vec+m*(p-1), tmp, tmp);
			part2 += dmat_norm1(p-1, tmp);
		}
        
		for (int l = 0; l<L; l++) {
			part3 += .5*a[l]*dmat_xAx(p-1, A+l*(p-1)*(p-1), X+l*(p-1));
		}
        part3 += dmat_dot((p-1)*L, X, b);
        
		return (part1+part2+part3);
	}
    
        // function to solve sub lasso problem: \sum_l (1/2 a_l x_l^T A_l x_l + b_l^T x_l + lambda1 * sum_j |x_lj|) + lambda2 * sum_{l l' in graph} |zeta_ll'| + mu / 2 || x_l - x_l' -zeta_ll' + tau_ll'/mu  ||^2
    void solve_lasso(double* X, double* tau_vec, double* zeta_vec, double* a, double* b, double* A, double* lambda1_vec, double* lambda2_vec, int *graph, double *mu, int p, int L, int NumEdges, double* EPS){
        int MAX_ITER = 200;
        double eps = EPS[0];
        vector<vector<int> > graph_adj(2*L,vector<int>(0)); // l-th entries: l1 m1, l2 m2 ....... li connected with l and is the m-th edges in graph.
        for (int m = 0; m < NumEdges; m++) {
            int l1 = graph[2*m];
            int l2 = graph[2*m+1];
            graph_adj[2*l1].push_back(l2);
            graph_adj[2*l1].push_back(m);
            graph_adj[2*l2+1].push_back(l1);
            graph_adj[2*l2+1].push_back(m);
        }
        double diff;
        
        vector<vector<int> > active_set_X(L,vector<int>(0));
        
        for (int iter = 0; iter < MAX_ITER; iter++) {
            
            double X_old[(p-1)*L];
            dmat_vcopy((p-1)*L, X, X_old);
            
                // circle once over X and zeta_vec //
                // update X //
            for (int l = 0; l < L; l++) {
                for (int k = 0; k < p-1; k++) {
                    double c_lk = 0.0;
                    double a_lk = 0.0;
                        // vector<int> tmp1 = graph_adj[2*l];
                        // vector<int> tmp2 = graph_adj[2*l+1];
                    int d_l = ((int)graph_adj[2*l].size()+(int)graph_adj[2*l+1].size()) / 2;
                    
                    for (int j = 0; j < p-1; j++) {
                        c_lk += X[l*(p-1)+j]*A[l*(p-1)*(p-1)+j*(p-1)+k];
                    }
                    c_lk = a[l]*(c_lk - X[l*(p-1)+k]*A[l*(p-1)*(p-1)+k*(p-1)+k]);
                    
                    if (graph_adj[2*l].size() > 0) {
                        for (int i = 0; i < graph_adj[2*l].size()/2;i++){
                                // int m = graph_adj[2*l][2*i+1];
                                // int l_prime = graph_adj[2*l][2*i];
                                // c_lk = c_lk - mu[0]*(X[l_prime*(p-1)+k] + zeta_vec[m*(p-1)+k]) + tau_vec[m*(p-1)+k];
                            c_lk = c_lk - mu[0]*(X[graph_adj[2*l][2*i]*(p-1)+k] + zeta_vec[graph_adj[2*l][2*i+1]*(p-1)+k]) + tau_vec[graph_adj[2*l][2*i+1]*(p-1)+k];
                        }
                    }
                    if (graph_adj[2*l+1].size() > 0) {
                        for (int i = 0; i < graph_adj[2*l+1].size()/2;i++){
                                // int m = graph_adj[2*l+1][2*i+1];
                                //int l_prime = graph_adj[2*l+1][2*i];
                                // c_lk = c_lk + mu[0]*(X[l_prime*(p-1)+k] + zeta_vec[m*(p-1)+k]) - tau_vec[m*(p-1)+k];
                            c_lk = c_lk - mu[0]*(X[graph_adj[2*l+1][2*i]*(p-1)+k] - zeta_vec[graph_adj[2*l+1][2*i+1]*(p-1)+k]) - tau_vec[graph_adj[2*l+1][2*i+1]*(p-1)+k];
                        }
                    }
                    a_lk = a[l]*A[l*(p-1)*(p-1)+k*(p-1)+k] + d_l * mu[0];
                    c_lk = - (c_lk + b[l*(p-1)+k]);
                    
                        // X[l*(p-1)+k] = sign(c_lk) * max(c_lk-lambda1_vec[l*(p-1)+k], 0) / a_lk;
                    soft_thred(1, X+l*(p-1)+k, &c_lk, &a_lk, lambda1_vec+l*(p-1)+k);
                    if (X[l*(p-1)+k] != 0) {
                        active_set_X[l].push_back(k);
                    }
                }
            }
                // update zeta //
            for (int m = 0; m < NumEdges; m++) {
                double c_m[p-1];
                    // int l1 = graph[2*m];
                    // int l2 = graph[2*m+1];
                    // dmat_waxpby(p-1, 1.0, X+l1*(p-1), -1.0, X+l2*(p-1), c_m);
                dmat_waxpby(p-1, mu[0], X+graph[2*m]*(p-1), -mu[0], X+graph[2*m+1]*(p-1), c_m);
                dmat_waxpby(p-1, 1.0, c_m, 1.0, tau_vec+m*(p-1), c_m);
                soft_thred(p-1, zeta_vec+m*(p-1), c_m, mu, lambda2_vec+m*(p-1));
            }
            
            /*
             for (int l = 0; l < L; l++) {
             for (int k = 0; k < (int)active_set_X[l].size(); k++) {
             printf("%d ", active_set_X[l][k]);
             }
             printf("\n");
             }
             
             print_dmatrix(X_old, 1, (p-1));
             print_dmatrix(X, 1, (p-1)); */
            
            if ((zero_agree((p-1)*L,X_old,X) == 1) && (iter > 0)) {
                    //printf("Lasso methods converges in %d iterations \n", iter);
                break;
            }
                // circle active set //
            for (int iter_active = 0; iter_active < MAX_ITER; iter_active++) {
                diff = 0.0;
                
                    // update X //
                for (int l = 0; l < L; l++) {
                    if (active_set_X[l].size()>0) {
                        for (int i = 0; i < active_set_X[l].size(); i++) {
                            int k = active_set_X[l][i];
                            double c_lk = 0.0;
                            double a_lk = 0.0;
                                // vector<int> tmp1 = graph_adj[2*l];
                                // vector<int> tmp2 = graph_adj[2*l+1];
                            int d_l = ((int)graph_adj[2*l].size()+(int)graph_adj[2*l+1].size()) / 2;
                            for (int jj = 0; jj < active_set_X[l].size(); jj++) {
                                c_lk += X[l*(p-1)+active_set_X[l][jj]]*A[l*(p-1)*(p-1)+active_set_X[l][jj]*(p-1)+k];
                            }
                            c_lk = a[l]*(c_lk - X[l*(p-1)+k]*A[l*(p-1)*(p-1)+k*(p-1)+k]);
                            
                            if(graph_adj[2*l].size() > 0){
                                for (int i = 0; i < graph_adj[2*l].size()/2;i++){
                                        // int m = graph_adj[2*l][2*i+1];
                                        // int l_prime = graph_adj[2*l][2*i];
                                        // c_lk = c_lk - mu[0]*(X[l_prime*(p-1)+k] + zeta_vec[m*(p-1)+k]) + tau_vec[m*(p-1)+k];
                                    c_lk = c_lk - mu[0]*(X[graph_adj[2*l][2*i]*(p-1)+k] + zeta_vec[graph_adj[2*l][2*i+1]*(p-1)+k]) + tau_vec[graph_adj[2*l][2*i+1]*(p-1)+k];
                                }
                            }
                            if (graph_adj[2*l+1].size() > 0) {
                                for (int i = 0; i < graph_adj[2*l+1].size()/2;i++){
                                        // int m = graph_adj[2*l+1][2*i+1];
                                        //int l_prime = graph_adj[2*l+1][2*i];
                                        // c_lk = c_lk + mu[0]*(X[l_prime*(p-1)+k] + zeta_vec[m*(p-1)+k]) - tau_vec[m*(p-1)+k];
                                    c_lk = c_lk - mu[0]*(X[graph_adj[2*l+1][2*i]*(p-1)+k] - zeta_vec[graph_adj[2*l+1][2*i+1]*(p-1)+k]) - tau_vec[graph_adj[2*l+1][2*i+1]*(p-1)+k];
                                }
                            }
                            a_lk = a[l]*A[l*(p-1)*(p-1)+k*(p-1)+k] + d_l * mu[0];
                            c_lk = - (c_lk + b[l*(p-1)+k]);
                            
                                // X[l*(p-1)+k] = sign(c_lk) * max(c_lk-lambda1_vec[l*(p-1)+k], 0) / a_lk;
                            double newx;
                            soft_thred(1, &newx, &c_lk, &a_lk, lambda1_vec+l*(p-1)+k);
                            double diff_new = newx - X[l*(p-1)+k];
                            diff = max(abs(diff_new),diff);
                            X[l*(p-1)+k] = newx;
                        }
                    }
                }
                
                    // update zeta //
                for (int m = 0; m < NumEdges; m++) {
                    double c_m[p-1];
                        // int l1 = graph[2*m];
                        // int l2 = graph[2*m+1];
                        // dmat_waxpby(p-1, 1.0, X+l1*(p-1), -1.0, X+l2*(p-1), c_m);
                    dmat_waxpby(p-1, mu[0], X+graph[2*m]*(p-1), -mu[0], X+graph[2*m+1]*(p-1), c_m);
                    dmat_waxpby(p-1, 1.0, c_m, 1.0, tau_vec+m*(p-1), c_m);
                    soft_thred(p-1, zeta_vec+m*(p-1), c_m, mu, lambda2_vec+m*(p-1));
                }
                
                if (diff < eps) {
                    break;
                }
            }
            
                // clear vectors objects on the heap //
            for (int l = 0; l < L; l++) {
                active_set_X[l].clear();
            }
        }
            // clear vectors in graph_adj //
        for (int j = 0; j < (int)graph_adj.size(); j++) {
            graph_adj[j].clear();
        }
        graph_adj.clear();
        active_set_X.clear();
        
    }
    
    
        // function to solve a sub-problem: \sum_l (1/2 a_l x_l^T A_l x_l + b_l^T x_l + lambda1 * sum_j |x_lj|) + lambda2 * sum_{l l' in graph} |x_l - x_l'|
    void solve_grouping(double *X, double *a,double *b,double *A,int *graph,double *lambda1_vec,double *lambda2_vec,int p,int L,int NumEdges,int col){
        
            // start Augmented lagrangian multiplier methods //
        double rho = 2.0;
        double mu = 3.0;
        double eps = 1e-6;
        double lasso_eps = 1e-6;
        int MAX_ITER = 40;
        double tau_vec[(p-1)*NumEdges];
        dmat_vset((p-1)*NumEdges, 0.0, tau_vec);
        double zeta_vec[(p-1)*NumEdges];
        for (int m = 0; m < NumEdges; m++) {
            dmat_waxpby(p-1, 1.0, X+(p-1)*graph[2*m], -1.0, X+(p-1)*graph[2*m+1], zeta_vec+(p-1)*m);
        }
        double prev = group_fun(X, a, b, A, graph, lambda1_vec, lambda2_vec, p, L, NumEdges);
        double cur, decre;
        for (int iter = 0; iter < MAX_ITER; iter++) {
            /*
             printf("\n new iteration starts: \n \n");
             print_dmatrix(X,1,(p-1)*L);
             print_dmatrix(tau_vec,1,(p-1)*NumEdges);
             print_dmatrix(zeta_vec,1,(p-1)*NumEdges);
             print_dmatrix(a, L, 1);
             print_dmatrix(b, (p-1)*L, 1);
             print_dmatrix(A, (p-1), (p-1)*L);
             print_dmatrix(lambda1_vec, (p-1)*L, 1);
             print_dmatrix(lambda2_vec, (p-1)*NumEdges, 1);
             print_imatrix(graph,2,NumEdges);
             */
            
                // update X and zeta_vec //
            solve_lasso(X,tau_vec,zeta_vec,a,b,A,lambda1_vec,lambda2_vec,graph,&mu,p,L,NumEdges,&lasso_eps);
            
                // convergence check starts ...//
            
                // feasibility check //
            double feas = 0.0;
            double tmp[p-1];
            for (int m = 0; m < NumEdges; m++) {
                dmat_waxpby(p-1, 1.0, X+graph[2*m]*(p-1), -1.0, X+graph[2*m+1]*(p-1), tmp);
                dmat_waxpby(p-1, 1.0, zeta_vec+m*(p-1), -1.0, tmp, tmp);
                feas += dmat_norm1(p-1, tmp);
            }
            feas = feas / max(dmat_norm1((p-1)*L, X),1);
            
                // function value check
            cur = group_fun(X, a, b, A, graph, lambda1_vec, lambda2_vec, p, L, NumEdges);
            decre = (prev - cur) / (abs(prev)+1.0);
            
            if (max(abs(decre),feas) < eps){
                    // printf("Augmented Lagrangian methods converges in %d iterations \n", iter);
                break;
            }
            prev = cur;
                // convergence check ends ...//
            
                // update the Lagrangian multipliers: tau_vec //
            for (int m = 0; m < NumEdges; m++) {
                dmat_waxpby(p-1, 1.0, tau_vec+(p-1)*m, mu, X+(p-1)*graph[2*m], tau_vec+(p-1)*m); //tau = tau + mu*xl
                dmat_waxpby(p-1, 1.0, tau_vec+(p-1)*m, -mu, X+(p-1)*graph[2*m+1], tau_vec+(p-1)*m); //tau = tau - mu*xl'
                dmat_waxpby(p-1, 1.0, tau_vec+(p-1)*m, -mu, zeta_vec+(p-1)*m, tau_vec+(p-1)*m); // tau = tau - mu*zeta_ll'
            }
                // update penalty: mu //
            mu = mu * rho;
        }
    }
    
		// solving convex version: min. \sum_{l = 1}^L n_l(-logdet(X_l) + tr(S_l * X_l)) + lambda1 * sum_{l=1}^L (n_l * \|X_l\|_{off}) + lambda2 * sum_{l ~ l'} \|X_l - X_l'\|_{off}
    void matrix_grouping_sub(double *covmat_inverse,double *covmat,double *S,double *lambda1,double *lambda2,double *tau,int *graph,double *nn,int p,int L,int NumEdges,int *lambda1_mat,int *lambda2_mat, double* eps_mat){
        double* A = (double*)malloc(sizeof(double)*(p-1)*(p-1)*L);
        double* covmat_inverse_old = (double*)malloc(sizeof(double)*p*p*L);
        int MAX_ITER = 100;
        for (int iter_col = 0; iter_col < MAX_ITER; iter_col++) {
            dmat_vcopy(p*p*L, covmat_inverse, covmat_inverse_old);
                // sweeping across columns of covmat_inverse //
            for (int col = 0; col < p; col++) {
                double a[L];
                for (int l = 0; l < L; l++) {
                    a[l] = 2*nn[l]*S[l * p * p + col * p + col];
                }
                double b[L * (p - 1)];
                for (int l = 0; l < L; l++){
                    for (int j = 0; j < p - 1; j++){
                        int j_index = ((j >= col)? (j+1):j);
                        b[l*(p-1)+j] = 2*nn[l]*S[l*p*p+col*p+j_index];
                    }
                }
                update_A(covmat, A, col, p, L);
                
                double X[(p-1)*L];
                for (int i = 0; i < p-1; i++) {
                    int i_index = ((i>=col)? (i+1):i);
                    for (int l = 0; l<L; l++) {
                        X[(p-1)*l + i] = covmat_inverse[l*p*p+col*p+i_index];
                    }
                }
                
                double lambda1_vec[(p-1)*L];
                double lambda2_vec[(p-1)*NumEdges];
                for (int j = 0; j < p-1; j++) {
                    int index_j = (j >= col)? (j+1):(j);
                    for (int l = 0; l < L; l++) {
                        // lambda1_vec[l*(p-1)+j] = lambda1[0]*lambda1_mat[l*p*p+col*p+index_j];
                        // different penalization for different sample sizes //
                        lambda1_vec[l*(p-1)+j] = lambda1[0] * nn[l] * lambda1_mat[l*p*p+col*p+index_j];
                    }
                    
                    for (int m = 0; m < NumEdges; m++) {
                        lambda2_vec[m*(p-1)+j] = lambda2[0] * lambda2_mat[m*p*p+col*p+index_j];
                    }
                }
                
                    // print_dmatrix(A, (p-1), (p-1)*L);
                    // print_dmatrix(X, (p-1)*L, 1);
                solve_grouping(X,a,b,A,graph,lambda1_vec,lambda2_vec,p,L,NumEdges,col);
                
                    // solve x_ljj -- the diagnol elements : 1/s_jjl + x_l-j^T (X_l)_{-j}^{-1} x_l-j
                for (int l = 0; l<L; l++) {
                    covmat_inverse[l*p*p+col*p+col] = 1/S[l * p * p + col * p + col] + dmat_xAx(p-1, A+l*(p-1)*(p-1), X+l*(p-1));
                }
                
                    // update the col-th column of every matrix in current_matrix
                for (int l = 0; l<L; l++) {
                    for (int i = 0; i<(p-1); i++) {
                        int i_index = ((i >= col)? (i+1):i);
                        covmat_inverse[l*p*p+col*p+i_index] = X[l*(p-1)+i];
                        covmat_inverse[l*p*p+i_index*p+col] = X[l*(p-1)+i];
                    }
                }
                
                update_cov(covmat,covmat_inverse, A, col, p, L);
            }
            
                // covmat_inverse_old = covmat_inverse - covmat_inverse_old //
            dmat_waxpby(p*p*L, 1.0, covmat_inverse, -1.0, covmat_inverse_old, covmat_inverse_old);
            double diff = dmat_norm1(p*p*L, covmat_inverse_old) / dmat_norm1(p*p*L, covmat_inverse);
                // printf("diff is %f \n ", diff);
                // printf("sweeping round %d... \n",iter_col);
            if (diff < eps_mat[0]) {
                    // printf("sweeping converges at %d iterations!! \n \n \n", iter_col);
                break;
            }
        }
        free(A); free(covmat_inverse_old);
    }
    
    void update_lambda(double *X,double *tau,int *graph,int p,int L,int NumEdges,int *lambda1_mat,int *lambda2_mat,int *converge){
            // update the l-th matrix with j-th col and k-th row (updates upper part)
        int num_diff = converge[0];
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < p; k++) {
                for (int l = 0; l < L; l++) {
                    if (abs(X[l*p*p+j*p+k]) > tau[0]) {
                        if (lambda1_mat[l*p*p+j*p+k]==1) {
                            lambda1_mat[l*p*p+j*p+k] = 0;
                            num_diff++;
                        }
                    } else {
                        if (lambda1_mat[l*p*p+j*p+k]==0) {
                            lambda1_mat[l*p*p+j*p+k] = 1;
                            num_diff++;
                        }
                    }
                }
                
                for (int m =0; m < NumEdges; m++) {
                    double tmp = X[graph[2*m]*p*p+j*p+k] - X[graph[2*m+1]*p*p+j*p+k];
                    tmp = abs(tmp);
                    if (tmp > tau[0]) {
                        if (lambda2_mat[m*p*p+j*p+k]==1) {
                            lambda2_mat[m*p*p+j*p+k] = 0;
                            num_diff++;
                        }
                    } else {
                        if (lambda2_mat[m*p*p+j*p+k]==0) {
                            lambda2_mat[m*p*p+j*p+k] = 1;
                            num_diff++;
                        }
                    }
                }
            }
        }
        converge[0] = (num_diff > 0)? 0:1;
    }
    
        // solving non-convex version: min. \sum_{l = 1}^L n_l(-logdet(X_l) + tr(S_l * X_l)) + lambda1 * sum_{l=1}^L (n_l * \|X_l\|_{off}) + lambda2 * sum_{l ~ l'} \|X_l - X_l'\|_{off}
    void matrix_grouping(double *covmat_inverse0, double *covmat_inverse_con0, double *covmat0, double *covmat_con0, double *S_bar,double *Lambda1,double *Lambda2, double *Tau,int *Graph,double *sample_size,int *pp,int *LL, int *NumOfEdges,double* eps_mat){
        int MAX_DC_ITER = 20;
        double *covmat_inverse = covmat_inverse0;
        double *covmat_inverse_con = covmat_inverse_con0;
        double *covmat = covmat0;
        double *covmat_con = covmat_con0;
        double *S = S_bar;
        double *lambda1 = Lambda1;
        double *lambda2 = Lambda2;
        double *tau = Tau;
        int *graph = Graph;
        double *nn = sample_size;
        int p = pp[0];
        int L = LL[0];
        int NumEdges = NumOfEdges[0];
        int *lambda1_mat = (int*) malloc(sizeof(int)*p*p*L);
        int *lambda2_mat = (int*) malloc(sizeof(int)*p*p*NumEdges);
        
            // set lambda value for convex penalties //
        dmat_iset(p*p*L, 1, lambda1_mat);
        dmat_iset(p*p*NumEdges,1,lambda2_mat);
        
        
            //print_imatrix(lambda1_mat,p,p*L);
            //print_imatrix(lambda2_mat,p,p*NumEdges);
        
            // get convex solutions //
        matrix_grouping_sub(covmat_inverse_con,covmat_con,S,lambda1,lambda2,tau,graph,nn,p,L,NumEdges,lambda1_mat,lambda2_mat,eps_mat);
        
            // copy convex solution as initial value for tlp solutions //
        dmat_vcopy(p*p*L, covmat_inverse_con, covmat_inverse);
        dmat_vcopy(p*p*L,covmat_con,covmat);
        
        int converge = 0;
        for (int dc_iter = 0; dc_iter < MAX_DC_ITER; dc_iter++){
                // update lambda1_mat and lambda2_mat for penalty //
            update_lambda(covmat_inverse, tau, graph, p, L, NumEdges, lambda1_mat, lambda2_mat, &converge);
            
                //print_imatrix(lambda1_mat,p,p*L);
                //print_imatrix(lambda2_mat,p,p*NumEdges);
            if (converge == 1){
                    // printf("DC converges in %d iterations \n", dc_iter);
                break;
            }
            
            matrix_grouping_sub(covmat_inverse,covmat,S,lambda1,lambda2,tau,graph,nn,p,L,NumEdges,lambda1_mat,lambda2_mat,eps_mat);
        }
        free(lambda1_mat); free(lambda2_mat);
        
    }
    
        // obtaining solution paths //
    void matrix_grouping_path(double *S_bar, double *covmat_inverse_path, double *covmat_path, 
		double *covmat_inverse_con_path, double *covmat_con_path,double *Lambda1_path,
			double *Lambda2_path, double *Tau,int* grid_lambda1, int* grid_lambda2, int *Graph,
				double *sample_size,int *pp,int *LL, int *NumOfEdges,double* eps_mat){
        int MAX_DC_ITER = 10;
        double *covmat_inverse = covmat_inverse_path;
        double *covmat_inverse_con = covmat_inverse_con_path;
        double *covmat = covmat_path;
        double *covmat_con = covmat_con_path;
        double *S = S_bar;
        double *Lambda1 = Lambda1_path;
        double *Lambda2 = Lambda2_path;
        int lambda1_grid = grid_lambda1[0];
        int lambda2_grid = grid_lambda2[0];
        double *tau = Tau;
        int *graph = Graph;
        double *nn = sample_size;
        int p = pp[0];
        int L = LL[0];
        int NumEdges = NumOfEdges[0];
        int *lambda1_mat = (int*) malloc(sizeof(int)*p*p*L);
        int *lambda2_mat = (int*) malloc(sizeof(int)*p*p*NumEdges);
        int index_diff = 0;
        for (int lambda1_id = 0; lambda1_id < lambda1_grid; lambda1_id++) {
            for (int lambda2_id = 0; lambda2_id < lambda2_grid ; lambda2_id++) {
                // printf("%f \t %f \n",Lambda1[lambda1_id], Lambda2[lambda2_id]);
                    // set lambda value for convex penalties //
                dmat_iset(p*p*L, 1, lambda1_mat);
                dmat_iset(p*p*NumEdges,1,lambda2_mat);
                index_diff = (lambda1_id*lambda2_grid + lambda2_id)*p*p*L;
                
                    // get convex solutions (initializing with previous solution)//
                if (index_diff != 0) {
                    dmat_vcopy(p*p*L,covmat_inverse_con+index_diff-p*p*L,covmat_inverse_con+index_diff);
                    dmat_vcopy(p*p*L,covmat_con+index_diff-p*p*L,covmat_con+index_diff);
                }
                matrix_grouping_sub(covmat_inverse_con+index_diff,covmat_con+index_diff,S,Lambda1+lambda1_id,Lambda2+lambda2_id,tau,graph,nn,p,L,NumEdges,lambda1_mat,lambda2_mat,eps_mat);
                
                    // copy convex solution as initial value for tlp solutions //
                dmat_vcopy(p*p*L, covmat_inverse_con+index_diff, covmat_inverse+index_diff);
                dmat_vcopy(p*p*L,covmat_con+index_diff,covmat+index_diff);
                
                int converge = 0;
                for (int dc_iter = 0; dc_iter < MAX_DC_ITER; dc_iter++){
                        // update lambda1_mat and lambda2_mat for penalty //
                    update_lambda(covmat_inverse+index_diff, tau, graph, p, L, NumEdges, lambda1_mat, lambda2_mat, &converge);
                    
                    if (converge == 1){
						// printf("\n DC converges in %d iterations \n", dc_iter);
                        break;
                    }
                    matrix_grouping_sub(covmat_inverse+index_diff,covmat+index_diff,S,Lambda1+lambda1_id,Lambda2+lambda2_id,tau,graph,nn,p,L,NumEdges,lambda1_mat,lambda2_mat,eps_mat);
                }
            }
        }
        
        free(lambda1_mat); free(lambda2_mat);
        
    }
	
	/*
		// .Call() implementation // 
	SEXP matrix_grouping_path_R(SEXP S_bar, SEXP sample_size, SEXP Lambda1.vec, SEXP Lambda2.vec, 
			SEXP Graph, SEXP Tau, SEXP MAX_iter, SEXP eps_mat) 
	{
		double *S, *Lambda1, *Lambda2, *tau, *nn, *graph;
		int MAX_DC_ITER = 20, lambda1_grid, lambda2_grid, p, L, NumEdges;
		S = REAL(S_bar); nn = REAL(sample_size); 
		Lambda1 = REAL(Lambda1.vec); Lambda2 = REAL(Lambda2.vec);
		graph = REAL(Graph); tau = REAL(Tau);
		SEXP Rdim,results;
		
	}*/
}

