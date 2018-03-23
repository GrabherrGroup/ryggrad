
#define FORCE_DEBUG

#include "../../src/base/CommandLineParser.h"
#include "../../src/base/FileParser.h"
#include "../../src/base/SVector.h"
#include "../../src/visual/Whiteboard.h"
#include "../../src/visual/Color.h"

#include "../../src/visual/Axes.h"

#include <iostream>

class Read
{
public:
  Read() {}
  Read(int id, int start, int stop, int partner, int ori) {
    m_start = start;
    m_stop = stop;
    m_partner = partner;
    m_ori = ori;
    m_id = id;
  }
  
  int ID() const {return m_id;}
  int Start() const {return m_start;}
  int Stop() const {return m_stop;}
  int Partner() const {return m_partner;}
  int Ori() const {return m_ori;}

private:
  int m_id;
  int m_start;
  int m_stop;
  int m_partner;
  int m_ori;
  
};


class Pos
{
 public:
  Pos() {
    x1 = x2 = y1 = y2 = -1;
    partner = -1;
    id = -1;
  }
  
  double x1;
  double x2;
  double y1;
  double y2;
  int partner;
  int id;
};

int main( int argc, char** argv )
{
  //No call to RunTime() made in order to allow clean ctrl-C exit.
 
  commandArg<string> aStringI1("-i","Layout file");
  commandArg<string> oStringI1("-o","Postscript file");
  commandArg<int> iCmd("-index","contig index", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Displays a transcript/contig");

  P.registerArg(aStringI1);
  P.registerArg(oStringI1);
  P.registerArg(iCmd);

  P.parse();

  string in = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(oStringI1);
  int index = P.GetIntValueFor(iCmd);
  
  double x_offset = 20;
  double y_offset = 20;

  int i, j;
  ns_whiteboard::whiteboard board;

  FlatFileParser parser;
  
  parser.Open(in);
  

  double x_max = 0;
  double y_max = 0;


  double x = 0;
  double y = 0;

  int k = 0;

  svec<Read> reads;

  bool bDo = false;
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "</CONTIG>") {
      if (bDo)
	break;
    }
  
    if (parser.AsString(0) == "<CONTIG>") {
      if (k == index) 
	bDo = true;
      k++;
      continue;
    }
    if (bDo) {
      if (parser.AsInt(5) >= 0) {
	reads.push_back(Read(parser.AsInt(0), 
			     parser.AsInt(2), 
			     parser.AsInt(4), 
			     parser.AsInt(5), 
			     parser.AsInt(1)));
      }
    }
  }


  svec<int> used;
  svec<int> part;
  used.resize(reads.isize());
  part.resize(reads.isize());

  svec<Pos> pos;
  pos.resize(reads.isize());
  
  int nn = 0;
  int fw = 0;
  int rc = 0;
  int correct = 0;
  for (i=0; i<reads.isize(); i++) {
    for (j=0; j<reads.isize(); j++) {
      if (reads[i].ID() == reads[j].Partner()) {
	part[i] = 1;
	part[j] = 1;
	nn++;	
	if (reads[i].Ori() == 1)
	  fw++;
	if (reads[i].Ori() == -11)
	  rc++;
	if (reads[i].Ori() != reads[j].Ori())
	  correct++;
      }
    }
  }
  cout << "Total: " << reads.isize();
  cout << " partnered: " << nn << " forward: " << fw << " reverse: " << rc << " correct: " << correct << endl;

  int d = reads.isize();

  svec<int> partner;
  partner.resize(reads.isize(), -1);
  
  for (i=0; i<reads.isize(); i++) {
    for (j=i+1; j<reads.isize(); j++) {
      if (reads[j].ID() == reads[i].Partner()) {
	partner[i] = j;
      }
    }
  }

  double scale = 0.5;
  
  while (d > 0) {
   
    int last = -1;
    x = 0.;
    d = reads.isize();

    for (i=0; i<reads.isize(); i++) {
      if (used[i] > 0) {
	d--;
	continue;
      }
      // if (d == 2)
      //cout << i << " " << reads[i].Pos() << " " << last << endl;

      if (reads[i].Start() > last) {
	int p = partner[i];

	if (p == -1) {
	  used[i] = 1;
	  d--;
	  continue;
	}

	int disc = 50;
	
	last = reads[p].Stop() + 5;
	double x = scale * reads[i].Start();
	double x1 = scale * (reads[i].Stop()-disc);
	if (x1 > x_max)
	  x_max = x1;

	pos[i].x1 = x;
	pos[i].x2 = x1;
	pos[i].y1 = y;
	pos[i].y2 = y;
	
	
	double r, g, b;
	r = g = b = 0.1;
	if (part[i] > 0) {
	  if (reads[i].Ori() == 1) {
	    r = 0.7;
	    g = 0.5;
	    b = 0.;
	  }
	  if (reads[i].Ori() == -1) {
	    r = 0.;
	    g = 0.5;
	    b = 0.7;
	  }
	}

	double x2 = x1;

	board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x_offset + x, y_offset + y),
					    ns_whiteboard::xy_coords(x_offset + x1, y_offset + y),
					    1.4, color(r, g, b)));
	
	x = scale * reads[p].Start();
	x1 = scale * (reads[p].Stop()-disc);
	if (x1 > x_max)
	  x_max = x1;

	pos[p].x1 = x;
	pos[p].x2 = x1;
	pos[p].y1 = y;
	pos[p].y2 = y;

	if (part[p] > 0) {
	  if (reads[p].Ori() == 1) {
	    r = 0.7;
	    g = 0.5;
	    b = 0.;
	  }
	  if (reads[p].Ori() == -1) {
	    r = 0.;
	    g = 0.5;
	    b = 0.7;
	  }
	}

	board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x_offset + x, y_offset + y),
					    ns_whiteboard::xy_coords(x_offset + x1, y_offset + y),
					    1.4, color(r, g, b)));

	// Link!!
	board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x_offset + x2, y_offset + y),
					    ns_whiteboard::xy_coords(x_offset + x, y_offset + y),
					    0.2, color(.1, .1, .1)));


	used[i]++;
	used[p]++;
	d--;
      }
    }
    y += 3.;
    //cout << d << endl;
  } 

  /*
  for (i=0; i<reads.isize(); i++) {
    for (j=i+1; j<reads.isize(); j++) {
      //cout << reads[j].ID() << " - " << reads[i].Partner() << endl;
      if (reads[j].ID() == reads[i].Partner()) {
	
	//cout << "Found partner" << endl;
	board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x_offset + pos[i].x1, y_offset + pos[i].y1),
					    ns_whiteboard::xy_coords(x_offset + pos[j].x2, y_offset + pos[j].y1),
					    1., color(.1, .1, .1)));
      }
    }
    }*/
  
  y_max = y;
 
  ofstream out(o.c_str());
  //cout << "MAX: " <<  x_max + 4 * x_offset + space << "\t" <<  y_max + 2 * y_offset + space << endl;
  ns_whiteboard::ps_display display(out, x_max + x_offset, y_max + y_offset);
  board.DisplayOn(&display);
 

  return 0;
}
