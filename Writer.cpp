#include "Writer.h"
using namespace std;

void Writer::write(string fileName, Model& model)
{
	int vs = model.vertices.size();
	int sz = model.faces.size();
	ofstream fout(fileName);
	fout << "ply" << endl;
	fout << "format ascii 1.0" << endl;
	fout << "element vertex " << vs << endl;
	fout << "property float x" << endl;
	fout << "property float y" << endl;
	fout << "property float z" << endl;
	fout << "element face " << sz << endl;
	fout << "property list uchar int vertex_indices" << endl;
	fout << "property uint8 red" << endl;
	fout << "property uint8 green" << endl;
	fout << "property uint8 blue" << endl;
	fout << "end_header" << endl;
	for (int i = 0; i < vs; ++i)
	{
		for (int j = 0; j < 3; ++j)
			fout << model.vertices[i].p.x[j] << " ";
		fout << endl;
	}
	for (int i = 0; i < sz; ++i)
	{
		fout << "3 ";
		for (int j = 0; j < 3; ++j)
			fout << model.faces[i].vids[j] << " ";
		int label = model.faces[i].label;
		fout << 60 * (label % 4 + 1) << " " << 80 * ((label + 1) % 3 + 1) << " " << 50 * ((label + 2) % 5 + 1) << endl;
	}
	fout.close();
}