#include "Reader.h"
using namespace std;

void Reader::read(string fileName, Model& model)
{
	ifstream fin(fileName);
	char tmp[256];
	char ch;
	double x, y, z;
	int v1, v2, v3;
	int fcount = 0;
	//int count = 0;
	while (fin.getline(tmp, 255))
	{
		//count++;
		//cout << "count = " << count << endl;
		sscanf(tmp, "%c", &ch);
		if (ch == '#')
			continue;
		else if (ch == 'v')
		{
			sscanf(tmp, "%c%lf%lf%lf", &ch, &x, &y, &z);
			Vertex v(x, y, z);
			model.vertices.push_back(v);
		}
		else if (ch == 'f')
		{
			sscanf(tmp, "%c%d%d%d", &ch, &v1, &v2, &v3);
			v1--;
			v2--;
			v3--;
			Face f(model.vertices[v1], model.vertices[v2], model.vertices[v3]);
			f.vids[0] = v1; f.vids[1] = v2; f.vids[2] = v3;
			model.faces.push_back(f);
			Edge e12(v1, v2, fcount);
			model.edges.push_back(e12);
			Edge e23(v2, v3, fcount);
			model.edges.push_back(e23);
			Edge e31(v3, v1, fcount);
			model.edges.push_back(e31);
			fcount++;
		}
	}
}