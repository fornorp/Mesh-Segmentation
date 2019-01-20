#pragma once
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <queue>

using namespace std;

#define M_PI 3.14159265358979323846
#define M_EPS 1e-12
#define M_INF 1e12

#define TH_DIST_RATIO 0.20
#define TH_ANGLE_DIFF 0.30
#define TH_REP_DIST_RATIO 0.1

struct Point
{
	double x[3];
	Point() {}
	Point(double xx, double yy, double zz)
	{
		x[0] = xx; x[1] = yy; x[2] = zz;
	}
	Point(const Point& p)
	{
		x[0] = p.x[0]; x[1] = p.x[1]; x[2] = p.x[2];
	}
	Point operator + (const Point& p) const
	{
		Point ans;
		for (int i = 0; i < 3; ++i)
			ans.x[i] = x[i] + p.x[i];
		return ans;
	}
	Point operator - (const Point& p) const
	{
		Point ans;
		for (int i = 0; i < 3; ++i)
			ans.x[i] = x[i] - p.x[i];
		return ans;
	}
	Point operator * (double d)
	{
		Point ans;
		for (int i = 0; i < 3; ++i)
			ans.x[i] = x[i] * d;
		return ans;
	}
	Point operator / (double d)
	{
		Point ans;
		for (int i = 0; i < 3; ++i)
			ans.x[i] = x[i] / d;
		return ans;
	}
}; 

typedef Point Vec;

double dot(const Vec& v1, const Vec& v2);
Point cross(const Vec& v1, const Vec& v2);
void normalizeV(Vec& v);

struct Vertex
{
	Point p;
	Vertex() {}
	Vertex(double x, double y, double z) : p(x, y, z) {}
	Vertex(const Vertex& v)
	{
		p = v.p;
	}
};

struct Neighbor
{
	int vids[2];
	int fid;
	double da, dg, d, a;
	Neighbor() {}
	Neighbor(int v1, int v2, int f, double da_, double dg_, double a_)
	{
		vids[0] = v1; vids[1] = v2; fid = f; da = da_; dg = dg_; a = a_;
	}
	Neighbor(const Neighbor& n)
	{
		vids[0] = n.vids[0]; vids[1] = n.vids[1];
		fid = n.fid;
		da = n.da; dg = n.dg; d = n.d; a = n.a;
	}
	void avg(double avg_da, double avg_dg, double delta)
	{
		d = (1 - delta) * da / avg_da + delta * dg / avg_dg;
		if (d < 0)
		{
			cout << d << " " << da / avg_da << " " << dg / avg_dg << endl;
		}
	}
};

struct Face
{
	int vids[3];
	vector<Neighbor> neighs;
	Vec normal;
	Point c;
	int label;
	Face() { label = 0; }
	Face(Vertex& v1, Vertex& v2, Vertex& v3)
	{
		c = (v1.p + v2.p + v3.p) / 3.0;
		normal = cross(v2.p - v1.p, v3.p - v1.p);
		normalizeV(normal);
		//if (dot(normal, c) < 0) normal = normal * (-1.0);
		label = 0;
	}
	Face(const Face& f)
	{
		vids[0] = f.vids[0]; vids[1] = f.vids[1]; vids[2] = f.vids[2];
		neighs = f.neighs;
		normal = f.normal;
		c = f.c;
		label = f.label;
	}
};

struct Edge
{
	int vids[2];
	int fid;
	Edge() {}
	Edge(int v1, int v2, int f)
	{
		vids[0] = v1; vids[1] = v2; fid = f;
	}
};

bool isConvex(const Face& fa, const Face& fb);
double faceGeoD(const Face& fa, const Face& fb, const Point& p0, const Point& p1);
double faceAngD(const Face& fa, const Face& fb, bool inverse=false);
double faceAng(const Face& fa, const Face& fb, bool inverse=false);

static double avg_da;
static double avg_dg;

struct Model
{
	std::vector<Vertex> vertices;
	std::vector<Face> faces;
	std::vector<Edge> edges;
	void findNeighFace(double delta)
	{
		double sum_da = 0.0, sum_dg = 0.0;
		int count = 0;
		for(int i = 0; i < edges.size(); ++i)
			for (int j = i + 1; j < edges.size(); ++j)
				if (edges[i].vids[0] == edges[j].vids[1] &&
					edges[i].vids[1] == edges[j].vids[0])
				{
					count++;
					double a = faceAng(faces[edges[i].fid], faces[edges[j].fid]);
					double da = faceAngD(faces[edges[i].fid], faces[edges[j].fid]);
					double dg = faceGeoD(faces[edges[i].fid], faces[edges[j].fid],
						vertices[edges[i].vids[0]].p, vertices[edges[i].vids[1]].p);
					sum_da += da;
					sum_dg += dg;
					Neighbor n1(edges[i].vids[0], edges[i].vids[1], edges[j].fid, da, dg, a);
					Neighbor n2(edges[j].vids[0], edges[j].vids[1], edges[i].fid, da, dg, a);
					faces[edges[i].fid].neighs.push_back(n1);
					faces[edges[j].fid].neighs.push_back(n2);
				}
				else if (edges[i].vids[0] == edges[j].vids[0] &&
					edges[i].vids[1] == edges[j].vids[1])
				{
					count++;
					double a = faceAng(faces[edges[i].fid], faces[edges[j].fid], true);
					double da = faceAngD(faces[edges[i].fid], faces[edges[j].fid], true);
					double dg = faceGeoD(faces[edges[i].fid], faces[edges[j].fid],
						vertices[edges[i].vids[0]].p, vertices[edges[i].vids[1]].p);
					sum_da += da;
					sum_dg += dg;
					Neighbor n1(edges[i].vids[0], edges[i].vids[1], edges[j].fid, da, dg, a);
					Neighbor n2(edges[j].vids[0], edges[j].vids[1], edges[i].fid, da, dg, a);
					faces[edges[i].fid].neighs.push_back(n1);
					faces[edges[j].fid].neighs.push_back(n2);
				}
		avg_da = sum_da / count;
		avg_dg = sum_dg / count;
		cout << "sum_da: " << sum_da << endl;
		cout << "count : " << count << endl;
		cout << "avg_da: " << avg_da << endl;
		for(int i = 0; i < faces.size(); ++i)
			for (int j = 0; j < faces[i].neighs.size(); ++j)
			{
				faces[i].neighs[j].avg(avg_da, avg_dg, delta);
				if (faces[i].neighs[j].d < -M_EPS)
					cout << "negative " << i << " " << j << " " << faces[i].neighs[j].d << endl;
			}
	}
};

static int types;
static double globalMaxD;
static double globalAvgD;
static int globalSz;

struct FlowEdge
{
	int id;
	double cap;
	double flow;
};

struct FlowState
{
	int id;
	int front;
	double inc;
};

class Graph
{
private:
	vector<vector<double> > weight;
	vector<vector<int> > neigh;
	vector<int> p;
	vector<int> heap;
	int tail;
	int sz;
	int s;

	vector<vector<FlowEdge> > net;
	vector<int> mask;
	int src, dst;

public:
	Graph() {};
	void build(Model& model)
	{
		sz = model.faces.size();
		dist.resize(sz);
		weight.resize(sz);
		neigh.resize(sz);
		net.resize(sz + 2);
		for (int i = 0; i < sz; ++i)
		{
			dist[i].resize(sz);
			for(int j = 0; j < sz; ++j)
				dist[i][j] = M_INF;
			dist[i][i] = 0;
		}
		heap.resize(sz + 2);
		p.resize(sz);
		for(int i = 0; i < sz; ++i)
			for (int j = 0; j < model.faces[i].neighs.size(); ++j)
			{
				weight[i].push_back(model.faces[i].neighs[j].d);
				neigh[i].push_back(model.faces[i].neighs[j].fid);
				FlowEdge fe;
				fe.cap = 1.0 / (1 + model.faces[i].neighs[j].da / avg_da);
				fe.flow = 0;
				fe.id = model.faces[i].neighs[j].fid;
				net[i].push_back(fe);
			}
		src = sz;
		dst = sz + 1;
	}
	void heapDown(int x)
	{
		while (x * 2 <= tail)
		{
			int y = x * 2;
			if (y < tail && dist[s][heap[y]] > dist[s][heap[y + 1]]) y++;
			if (dist[s][heap[x]] > dist[s][heap[y]])
			{
				int tmp = heap[x];
				heap[x] = heap[y];
				heap[y] = tmp;
				p[heap[x]] = x;
				p[heap[y]] = y;
				x = y;
			}
			else break;
		}
	}
	void heapUp(int x)
	{
		while (x > 1)
		{
			int y = x / 2;
			if (dist[s][heap[x]] < dist[s][heap[y]])
			{
				int tmp = heap[x];
				heap[x] = heap[y];
				heap[y] = tmp;
				p[heap[x]] = x;
				p[heap[y]] = y;
				x = y;
			}
			else break;
		}
	}
	void solve()
	{
		/*int head = 0, tail = 0;
		vector<int> q;
		vector<int> fill(sz);
		q.push_back(2074);
		fill[2074] = 1;
		while (head <= tail)
		{
			int now = q[head++];
			for (int i = 0; i < neigh[now].size(); ++i)
			{
				int next = neigh[now][i];
				if (!fill[next])
				{
					tail++;
					q.push_back(next);
					fill[next] = 1;
				}
			}
		}
		cout << "flood fill " << tail + 1 << " points" << endl;*/
		for (s = 0; s < sz; ++s)
		{
			int count = 0;
			for (int i = 0; i < sz; ++i)
				p[i] = 0;
			tail = 1;
			heap[tail] = s;
			p[s] = 1;
			while (tail > 0)
			{
				int now = heap[1];
				heap[1] = heap[tail];
				tail--;
				count++;
				heapDown(1);
				for (int i = 0; i < neigh[now].size(); ++i)
				{
					int next = neigh[now][i];
					double d = weight[now][i];
					if (dist[s][next] > dist[s][now] + d)
					{
						dist[s][next] = dist[s][now] + d;
						if (p[next] == 0) p[next] = ++tail;
						heap[p[next]] = next;
						heapUp(p[next]);
						heapDown(p[next]);
					}
				}
			}
		}
	}
	void initFlow(vector<int>& m)
	{
		mask = m;
		int netSz = 0;
		net[src].clear();
		for (int i = 0; i < sz; ++i)
			if (mask[i] == 1)
			{
				FlowEdge fe;
				fe.cap = M_INF;
				fe.flow = 0.0;
				fe.id = i;
				net[src].push_back(fe);
				netSz++;
			}
			else if (mask[i] == 2)
			{
				int last = net[i].size() - 1;
				if (net[i][last].id != dst)
				{
					FlowEdge fe;
					fe.cap = M_INF;
					fe.flow = 0.0;
					fe.id = dst;
					net[i].push_back(fe);
				}
				netSz++;
			}
			else
			{
				if (mask[i] == 3) netSz++;
				int last = net[i].size() - 1;
				if (net[i][last].id == dst)
					net[i].pop_back();
			}
		cout << "net size: " << netSz << endl;
		mask.push_back(4);
		mask.push_back(4);
		for (int i = 0; i < sz + 2; ++i)
			if (mask[i])
				for (int j = 0; j < net[i].size(); ++j)
					net[i][j].flow = 0.0;
	}
	void flow(vector<int>& m)
	{
		vector<FlowState> q;
		vector<bool> visit(sz + 2);
		int head = 0, tail = 0;
		int steps = 0;
		double sumFlow = 0;
		while (true)
		{
			steps++;
			head = tail = 0;
			q.clear();
			for (int i = 0; i < visit.size(); i++)
				visit[i] = false;
			FlowState fs;
			fs.id = src;
			fs.front = 0;
			fs.inc = M_INF;
			q.push_back(fs);
			visit[src] = true;
			while (head <= tail && !visit[dst])
			{
				FlowState now = q[head];
				for (int i = 0; i < net[now.id].size(); ++i)
				{
					int next = net[now.id][i].id;
					if (!mask[next]) continue;
					if (visit[next]) continue;
					double inc = net[now.id][i].cap - net[now.id][i].flow;
					if (inc < M_EPS) continue;
					fs.id = next;
					fs.front = head;
					if (now.inc > inc) fs.inc = inc;
					else fs.inc = now.inc;
					tail++;
					q.push_back(fs);
					visit[next] = true;
					if (next == dst) break;
				}
				head++;
			}
			if (!visit[dst])
			{
				for (int i = 1; i <= tail; ++i)
					mask[q[i].id] = 1;
				cout << "tail = " << tail << endl;
				for (int i = 0; i < sz; ++i)
					if (mask[i] == 3) mask[i] = 2;
				break;
			}
			double inc = q[tail].inc;
			sumFlow += inc;
			while (tail != 0)
			{
				int to = q[tail].id;
				int from = q[q[tail].front].id;
				if (inc > 1e8) cout << from << "->" << to << endl;
				for (int i = 0; i < net[from].size(); ++i)
					if (net[from][i].id == to)
						net[from][i].flow += inc;
				for (int i = 0; i < net[to].size(); ++i)
					if (net[to][i].id == from)
						net[to][i].flow -= inc;
				tail = q[tail].front;
			}
		}
		cout << "steps: " << steps << ", total flow: " << sumFlow << endl;
		/*int node = 2205;
		cout << "color: " << mask[node] << endl;
		cout << visit[node] << endl;
		for (int i = 0; i < net[node].size(); ++i)
			cout << net[node][i].cap << " " << net[node][i].flow << " " << net[node][i].id << endl;
		for (int i = 0; i < sz; ++i)
		{
			if (!mask[i]) continue;
			for (int j = 0; j < net[i].size(); ++j)
				if (net[i][j].id == node)
				{
					cout << i << " " << net[i][j].cap << " " << net[i][j].flow << endl;
				}
		}*/
		m = mask;
	}
	vector<vector<double> > dist;
};

class Solver
{
public:
	Solver() {}
	Solver(Model* m, Graph* g, double e = 0.0) { pModel = m; pGraph = g; eps = e; }
	void computeProb()
	{
		for (int i = 0; i < sz; ++i)
		{
			bool isRep = false;
			for(int j = 0; j < num; ++j)
				if (used[i] == reps[j])
				{
					isRep = true;
					for (int k = 0; k < num; ++k)
						prob[k][i] = 0;
					prob[j][i] = 1;
					break;
				}
			if (isRep) continue;
			double sum = 0;
			for (int j = 0; j < num; ++j)
				sum += 1.0 / pGraph->dist[used[i]][reps[j]];
			for (int j = 0; j < num; ++j)
				prob[j][i] = 1.0 / pGraph->dist[used[i]][reps[j]] / sum;
		}
	}
	void recomputeProb()
	{
		vector<int> typeNum(num);
		vector<vector<double> > distNew(num);
		for (int i = 0; i < num; ++i)
		{
			typeNum[i] = 0;
			distNew[i].resize(sz);
			for (int j = 0; j < sz; ++j)
				distNew[i][j] = 0;
		}
		for (int i = 0; i < sz; ++i)
		{
			int tmp = pModel->faces[used[i]].label - types;
			if (tmp < num) typeNum[tmp]++;
			else continue;
			for (int j = 0; j < sz; ++j)
			{
				distNew[tmp][j] += pGraph->dist[used[i]][used[j]];
			}
		}
		for (int i = 0; i < num; ++i)
		{
			for (int j = 0; j < sz; ++j)
				if (typeNum[i] > 0) distNew[i][j] /= typeNum[i];
				else distNew[i][j] = M_INF;
		}
		for (int i = 0; i < sz; ++i)
		{
			double sum = 0;
			for (int j = 0; j < num; ++j)
				sum += 1.0 / (distNew[j][i] + M_EPS);
			for (int j = 0; j < num; ++j)
				prob[j][i] = 1.0 / (distNew[j][i] + M_EPS) / sum;
		}
		//exit(1);
	}
	bool updateRep()
	{
		bool change = false;
		for (int i = 0; i < num; ++i)
		{
			double minSum = 0.0;
			for (int j = 0; j < sz; ++j)
				minSum += prob[i][j] * pGraph->dist[reps[i]][used[j]];
			for (int j = 0; j < sz; ++j)
			{
				//if (prob[i][j] < 0.5 / num) continue;
				double sum = 0;
				for (int k = 0; k < sz; ++k)
					sum += prob[i][k] * pGraph->dist[used[j]][used[k]];
				if (sum < minSum - M_EPS)
				{
					if (reps[i] != used[j]) change = true;
					minSum = sum;
					reps[i] = used[j];
				}
			}
		}
		return change;
	}
	void process()
	{
		int maxIter = 20;
		mask.resize(num);
		for (int i = 0; i < maxIter; ++i)
		{
			cout << "iter = " << i << endl;
			for (int j = 0; j < num; ++j)
				cout << reps[j] << " ";
			cout << endl;
			computeProb();
			assign();
			recomputeProb();
			if (!updateRep()) break;
		}
		cout << "converge" << endl;
		assign();
		cout << "after assign 0" << endl;
		recomputeProb();
		cout << "after recompute" << endl;
		assign();
		cout << "after assign 1" << endl;
		cut();
		types += num;
		this->expand();
	}
	void assign()
	{
		for (int i = 0; i < num; ++i)
		{
			int repeat = 0;
			for (int j = 0; j < i; ++j)
				if (reps[j] == reps[i])
					repeat = 1;
			mask[i] = 1-repeat;
		}
		for (int j = 0; j < num; ++j)
			cout << mask[j] << " ";
		cout << endl;
		int c[20] = {0};
		//double diff = eps / num;
		//double diff = eps / num * globalAvgD / localAvgD;
		double diff = eps * globalAvgD / localAvgD;
		if (diff > 0.2) diff = 0.2;
		diff = eps;
		cout << "diff = " << diff << endl;
		for (int i = 0; i < sz; i++)
		{
			double maxA = -1.0, maxB = -1.0;
			int labelA = 0, labelB = 0;
			for (int j = 0; j < num; j++)
			{
				if (!mask[j]) continue;
				if (maxA < prob[j][i])
				{
					maxB = maxA;
					labelB = labelA;
					maxA = prob[j][i];
					labelA = j;
				}
				else if (maxB < prob[j][i])
				{
					maxB = prob[j][i];
					labelB = j;
				}
			}
			if (maxA - maxB > diff)
			{
				pModel->faces[used[i]].label = types+labelA;
				c[labelA]++;
			}
			else
			{
				//pModel->faces[used[i]].label = (labelA+2)*(labelB+2);
				pModel->faces[used[i]].label = 1024 + (labelA*num) + labelB;
			}
		}
		cout << "split: " << sz << "->";
		for(int j = 0; j < num; ++j)
			cout << c[j] << " ";
		cout << endl;
		//types += num;
	}
	void cut()
	{
		vector<int> mm(globalSz);
		vector<int> mmNew;
		for (int i = 0; i < num; ++i)
		{
			if (!mask[i]) continue;
			for (int j = i + 1; j < num; ++j)
			{
				if (!mask[j]) continue;
				for (int k = 0; k < globalSz; ++k) mm[k] = 0;
				int nowLabelA = 1024 + i*num + j;
				int nowLabelB = 1024 + j*num + i;
				int fuzzy = 0;
				for(int k = 0; k < sz; ++k)
					if (pModel->faces[used[k]].label == nowLabelA || pModel->faces[used[k]].label == nowLabelB)
					{
						fuzzy++;
						mm[used[k]] = 3;
						for (int x = 0; x < pModel->faces[used[k]].neighs.size(); ++x)
						{
							int y = pModel->faces[used[k]].neighs[x].fid;
							if (pModel->faces[y].label == types + i) mm[y] = 1;
							else if (pModel->faces[y].label == types + j) mm[y] = 2;
						}
					}
				cout << "fuzzy: " << fuzzy << endl;
				pGraph->initFlow(mm);
				cout << "after init flow" << endl;
				pGraph->flow(mmNew);
				cout << "after flow" << endl;
				for (int k = 0; k < sz; ++k)
					if (mmNew[used[k]] == 1) pModel->faces[used[k]].label = types + i;
					else if (mmNew[used[k]] == 2) pModel->faces[used[k]].label = types + j;
					else if (mmNew[used[k]] != 0) cout << "wrong " << mmNew[used[k]] << endl;

				for (int k = 0; k < sz; ++k)
				{
					bool checked = false;
					for (int x = 0; x < pModel->faces[used[k]].neighs.size(); ++x)
					{
						int y = pModel->faces[used[k]].neighs[x].fid;
						if (pModel->faces[used[k]].label == pModel->faces[y].label)
							checked = true;
					}
					if (!checked)
						cout << "isolated face: " << used[k] << endl;
				}
			}
		}
	}
	double avgDistance()
	{
		double sum = 0.0;
		for (int i = 0; i < sz; ++i)
			for (int j = i + 1; j < sz; ++j)
				sum += pGraph->dist[used[i]][used[j]];
		sum /= sz*(sz - 1) / 2;
		return sum;
	}
	double angDifference()
	{
		double maxA = 0.0, minA = M_PI;
		for (int i = 0; i < sz; ++i)
		{
			int iLabel = pModel->faces[used[i]].label;
			for (int j = 0; j < pModel->faces[used[i]].neighs.size(); ++j)
			{
				int k = pModel->faces[used[i]].neighs[j].fid;
				if (k < used[i]) continue;
				int kLabel = pModel->faces[k].label;
				if (iLabel == kLabel)
				{
					double ang = pModel->faces[used[i]].neighs[j].a;
					if (ang < minA) minA = ang;
					if (ang > maxA) maxA = ang;
				}
			}
		}
		cout << maxA << " " << minA << endl;
		return maxA - minA;
	}
	double maxDistance()
	{
		double maxD = 0.0;
		for (int i = 0; i < sz; ++i)
			for (int j = i + 1; j < sz; ++j)
				if (maxD < pGraph->dist[used[i]][used[j]])
				{
					if (num == 2) 
					{
						reps[0] = used[i];
						reps[1] = used[j];
					}
					maxD = pGraph->dist[used[i]][used[j]];
				}
		cout << "in maxDistance " << maxD << endl;
		if (level == 0)
		{
			globalMaxD = maxD;
		}
		return maxD;
	}
	double maxPatchDistance()
	{
		double maxD = 0;
		for(int i = 0; i < num; ++i)
			for (int j = i + 1; j < num; ++j)
			{
				if (!mask[j]) continue;
				if (pGraph->dist[reps[i]][reps[j]] > maxD)
					maxD = pGraph->dist[reps[i]][reps[j]];
			}
		cout << "maxD = " << maxD << endl;
		cout << "global " << globalMaxD << endl;
		return maxD;
	}
	virtual void expand() { cout << "!!!" << endl; }
	virtual void init(vector<int>& u, int l) { cout << "???" << endl; }
	Model* pModel;
	Graph* pGraph;
	double eps;
	vector<int> used;
	vector<int> reps;
	vector<vector<double> > prob;
	vector<Solver*> sons;
	vector<bool> mask;
	int num;
	int sz;
	int level;
	double localAvgD;
	double localMaxD;
	double angRange;
};

class SolverTwo : public Solver
{
public:
	SolverTwo() {}
	SolverTwo(Model* m, Graph* g, double e = 0.0) : Solver(m, g, e) {}
	void init(vector<int>& u, int l)
	{
		level = l;
		used = u;
		num = 2;
		reps.resize(num);
		prob.resize(num);
		sz = used.size();
		for (int i = 0; i < num; ++i)
			prob[i].resize(sz);
		localMaxD = maxDistance();
		if (level == 0)
			globalAvgD = avgDistance();
		localAvgD = avgDistance();
		angRange = angDifference();
	}
	void expand()
	{
		double maxD = maxPatchDistance();
		cout << "minPatchD = " << maxD / globalMaxD << endl;
		if (level > 2 || maxD / globalMaxD < TH_REP_DIST_RATIO) return;
		for (int i = 0; i < num; ++i)
		{
			Solver* s = new SolverTwo(pModel, pGraph, eps);
			vector<int> uNew;
			for (int j = 0; j < sz; ++j)
				if ((pModel->faces[used[j]].label)%num == i)
					uNew.push_back(used[j]);
			s->init(uNew, level + 1);
			sons.push_back(s);
		}
		for (int i = 0; i < num; ++i)
		{
			cout << "maxD = " << sons[i]->localMaxD / globalMaxD 
				 << " angD = " << sons[i]->angRange << endl;
			if (sons[i]->localMaxD / globalMaxD > TH_DIST_RATIO && 
				sons[i]->angRange > TH_ANGLE_DIFF)
				sons[i]->process();
		}
	}
};

class SolverK : public Solver
{
public:
	SolverK() {}
	SolverK(Model* m, Graph* g, double e = 0.0) : Solver(m, g, e) {}
	void init(vector<int>& u, int l)
	{
		level = l;
		used = u;
		sz = used.size();
		double minSumD = M_INF;
		int choice = 0;
		for (int i = 0; i < sz; ++i)
		{
			double sumD = 0;
			for (int j = 0; j < sz; ++j)
				sumD += pGraph->dist[used[i]][used[j]];
			if (sumD < minSumD)
			{
				minSumD = sumD;
				choice = used[i];
			}
		}
		reps.push_back(choice);
		vector<double> G;
		double maxD = 0;
		for (int i = 1; i < 20; ++i)
		{
			maxD = 0;
			choice = 0;
			for (int j = 0; j < sz; ++j)
			{
				double minD = M_INF;
				for (int k = 0; k < i; ++k)
					if (minD > pGraph->dist[reps[k]][used[j]])
						minD = pGraph->dist[reps[k]][used[j]];
				if (minD > maxD)
				{
					maxD = minD;
					choice = used[j];
				}
			}
			reps.push_back(choice);
			G.push_back(maxD);
		}
		maxD = G[0] - G[1];
		num = 2;
		for (int i = 0; i < G.size(); ++i)
			cout << G[i] << " ";
		cout << endl;
		for(int i = 1; i < 18; ++i)
			if (maxD < G[i] - G[i + 1])
			{
				maxD = G[i] - G[i + 1];
				num = i + 2;
			}
		cout << "num = " << num << endl;
		//if (level == 0) num = 10;
		//else num = 2;
		while (reps.size() > num) reps.pop_back();
		prob.resize(num);
		for (int i = 0; i < num; ++i)
			prob[i].resize(sz);
		localMaxD = maxDistance();
		if (level == 0)
			globalAvgD = avgDistance();
		localAvgD = avgDistance();
		angRange = angDifference();
	}
	void expand()
	{
		double maxD = maxPatchDistance();
		cout << "minPatchD = " << maxD / globalMaxD << endl;
		if (level > 0 || maxD / globalMaxD < TH_REP_DIST_RATIO) return;
		for (int i = 0; i < num; ++i)
		{
			Solver* s = NULL;
			if (mask[i])
			{
				s = new SolverK(pModel, pGraph, eps);
				vector<int> uNew;
				for (int j = 0; j < sz; ++j)
					if ((pModel->faces[used[j]].label) % num == i)
						uNew.push_back(used[j]);
				s->init(uNew, level + 1);
			}
			sons.push_back(s);
		}
		for (int i = 0; i < num; ++i)
		{
			if (!mask[i]) continue;
			cout << "avgD = " << sons[i]->localAvgD / globalAvgD
				<< " angD = " << sons[i]->angRange << endl;
			if (sons[i]->localAvgD / globalAvgD > TH_DIST_RATIO &&
				sons[i]->angRange > TH_ANGLE_DIFF)
				sons[i]->process();
		}
	}
};