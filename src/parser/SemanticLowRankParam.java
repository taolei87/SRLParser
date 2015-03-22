package parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import utils.SVD;
import utils.Utils;

public class SemanticLowRankParam {
	public int N, M, D, L;
	public ArrayList<MatrixEntry> list;
	
	public SemanticLowRankParam(Parameters parameters) 
	{
		N = parameters.N2;
		M = parameters.M2;
		D = parameters.D2;
		L = parameters.L2;
		list = new ArrayList<MatrixEntry>();
	}

	
	
	public void add(int x, int y, int z, int r, double v)
	{
		list.add(new MatrixEntry(x, y, z, r, v));
	}
	
	public void decompose(Parameters params)
	{
		int maxRank = params.U2.length;
		
		int MAXITER=1000;
		double eps = 1e-6;
		Random rnd = new Random(0);
		for (int i = 0; i < maxRank; ++i) {
			double[] u = new double[N], v = new double[M], w = new double[D], r = new double[L];
			for (int j = 0; j < M; ++j)
				v[j] = rnd.nextDouble() - 0.5;
			for (int j = 0; j < D; ++j)
				w[j] = rnd.nextDouble() - 0.5;
			for (int j = 0; j < L; ++j)
				r[j] = rnd.nextDouble() - 0.5;
			Utils.normalize(v);
			Utils.normalize(w);
			Utils.normalize(r);
			int iter = 0;
			double norm = 0.0, lastnorm = Double.POSITIVE_INFINITY;
			for (iter = 0; iter < MAXITER; ++iter) {
				
				// u = <T,-,v,w,r>
				for (int j = 0; j < N; ++j)
					u[j] = 0;
				for (MatrixEntry e : list) {
					u[e.x] += e.value * v[e.y] * w[e.z] * r[e.r];
				}
				for (int j = 0; j < i; ++j) {
					double dot = Utils.dot(v, params.V2[j]) 
							   * Utils.dot(w, params.W2[j])
							   * Utils.dot(r, params.X2[j]);
					for (int k = 0; k < N; ++k)
						u[k] -= dot*params.U2[j][k];
				}
				Utils.normalize(u);
				
				// v = <T,u,-,w,r>
				for (int j = 0; j < M; ++j)
					v[j] = 0;
				for (MatrixEntry e : list) {
					v[e.y] += e.value * u[e.x] * w[e.z] * r[e.r];
				}
				for (int j = 0; j < i; ++j) {
					double dot = Utils.dot(u, params.U2[j]) 
							   * Utils.dot(w, params.W2[j])
							   * Utils.dot(r, params.X2[j]);
					for (int k = 0; k < M; ++k)
						v[k] -= dot*params.V2[j][k];
				}
				Utils.normalize(v);
				
				// w = <T,u,v,-,r>
				for (int j = 0; j < D; ++j)
					w[j] = 0;
				for (MatrixEntry e : list) {
					w[e.z] += e.value * u[e.x] * v[e.y] * r[e.r];
				}
				for (int j = 0; j < i; ++j) {
					double dot = Utils.dot(u, params.U2[j]) 
							   * Utils.dot(v, params.V2[j])
							   * Utils.dot(r, params.X2[j]);
					for (int k = 0; k < D; ++k)
						w[k] -= dot*params.W2[j][k];
				}
				Utils.normalize(w);
				
				// r = <T,u,v,w,->
				for (int j = 0; j < L; ++j)
					r[j] = 0;
				for (MatrixEntry e : list) {
					r[e.r] += e.value * u[e.x] * v[e.y] * w[e.z];
				}
				for (int j = 0; j < i; ++j) {
					double dot = Utils.dot(u, params.U2[j]) 
							   * Utils.dot(v, params.V2[j])
							   * Utils.dot(w, params.W2[j]);
					for (int k = 0; k < L; ++k)
						r[k] -= dot*params.X2[j][k];
				}			
				norm = Math.sqrt(Utils.squaredSum(r));
				if (lastnorm != Double.POSITIVE_INFINITY && Math.abs(norm-lastnorm) < eps)
					break;
				lastnorm = norm;
			}
			if (iter >= MAXITER) {
				System.out.printf("\tWARNING: Power method didn't converge." +
						"R=%d sigma=%f%n", i, norm);
			}
			if (Math.abs(norm) <= eps) {
				System.out.printf("\tWARNING: Power method has nearly-zero sigma. R=%d%n",i);
			}
			System.out.printf("\t%.2f", norm);
			params.U2[i] = u;
			params.V2[i] = v;
			params.W2[i] = w;
			params.X2[i] = r;
		}
		
		for (int i = 0; i < maxRank; ++i) {
			params.totalU2[i] = params.U2[i].clone();
			params.totalV2[i] = params.V2[i].clone();
			params.totalW2[i] = params.W2[i].clone();
			params.totalX2[i] = params.X2[i].clone();
		}
	}
}

class MatrixEntry
{
	public int x, y, z, r;
	public double value;
	
	public MatrixEntry(int _x, int _y, int _z, int _r, double _value)
	{
		x = _x;
		y = _y;
		z = _z;
		r = _r;
		value = _value;
	}
}

class MatrixEntryComparator implements Comparator<MatrixEntry>
{

	@Override
	public int compare(MatrixEntry o1, MatrixEntry o2) {
		if (o1.y != o2.y)
			return o1.y - o2.y;
		else
			return o1.x - o2.x;
	}
	
}