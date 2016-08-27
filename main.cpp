#include <stdio.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <list>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

//#include <ANN/ANN.h>   // ANN declarations

//#define DEBUG
//#define DEBUG_ALFA
//#define DEBUG_SHOW_S
//#define DEBUG_RBNN

using namespace std;


//------ ANN --------
int     k = 10;      // number of nearest neighbors
int     dim = 2;	// dimension
double  eps = 0;	// error bound
int	maxPts = 1000;	// maximum number of data points
//------ ANN --------

struct Point {
	double x;
	double y;
    double r;   // ray of scan
    bool hasCluster;
    int clusterIndex;
	int index;
};

const double PI = 3.14159265358979;

int scan_count=0;

FILE *f;

GLUquadricObj *diskQuadric;

template <class T>
T mabs ( T a )
{
  return sqrt(a*a);
}

float transZ = -10.0;

unsigned char segment = 0;

typedef struct _RGB_ {
  GLfloat r,g,b;
} RGBf;

RGBf seg_color[11];

typedef struct _seg_bound_ {
  int begin,end;
  vector<Point> p;
} seg_bound;

vector<seg_bound> segment_boundary;

float odo_x,odo_y,odo_phi;

bool polar_view = false;

int ALFA_line_count;

int th2 = 5; //num. point of line

double th1 = 0.1, th1_min = 0.2;
double m_RMSEth=0.1;

double Sx[700],Sy[700],Sxy[700],Sxx[700],Syy[700];

void BoundarySearch(seg_bound seg)
{
  cout << "\n\n---------- Boundary Search ----------\n\n" << endl;
  cout << "Seg Start: " << 0 << "\tSeg End: " << seg.p.size() << endl;  
  
  int k_min=0, k_max=seg.p.size();
  
  while ( (k_max - k_min) > th2 )
  {
    int b1 = k_min; int b2 = k_max;
    int a_ = round((b1+b2)/2.0);
    
#ifdef DEBUG_ALFA    
    cout << "a_: " << a_ << "  b1: " << b1 << "  b2: "<< b2 << endl;
#endif

    double RMSE_old = th1_min; th1 = 8.0;
    
    
    while ( (b2 - b1) > 0 )
    {
      int n=a_-k_min+1;
      
      // calculate LEAST SQUARE LINE FITTING CRITERIA      
      double SST = (Syy[a_+1]-Syy[k_min])-(Sy[a_+1]-Sy[k_min])*(Sy[a_+1]-Sy[k_min])/n;
      double SST_= (Sxx[a_+1]-Sxx[k_min])-(Sx[a_+1]-Sx[k_min])*(Sx[a_+1]-Sx[k_min])/n;
      double TMP = (n*(Sxy[a_+1]-Sxy[k_min])-(Sx[a_+1]-Sx[k_min])*(Sy[a_+1]-Sy[k_min]));
      double R2  = (TMP*TMP)/((n*SST_)*(n*SST));
      
#ifdef DEBUG_ALFA
      cout << "----------------SST: "  << SST << endl;
      cout << "----------------SST_: " << SST_ << endl;
      cout << "----------------TMP: "  << TMP << endl;
      cout << "----------------R2: "   << R2 << endl;
#endif

      double SSE,MSE,RMSE,m_RMSE;
      SSE=MSE=RMSE=m_RMSE=0.0;
      
      if (SST>SST_)
        SSE = (1-R2)*SST_;
      else
        SSE = (1-R2)*SST;

      MSE = SSE/n;
      RMSE = sqrt(MSE);
      
#ifdef DEBUG_ALFA
      cout << "RMSE: " << RMSE << "\tMSE: " << MSE << "\tSSE: " << SSE << "\tn: " << n << endl;
#endif

      // calculate slop of RMSE
      m_RMSE = (RMSE - max(th1_min,RMSE_old))/(a_ - b1);

      // check goodness of fit
      int bad_fit;
      bad_fit = (RMSE > th1) || (mabs(m_RMSE) > m_RMSEth);
      
      if  (bad_fit)
      {
      
#ifdef DEBUG_ALFA      
        cout << "[[bad fit RMSE: " << RMSE << " m_RMSE: " << mabs(m_RMSE);
#endif
        // cut a part of points for checking in next step
        b2 = a_-1;
        a_ = floor((b1 + b2)/2.0);
	
#ifdef DEBUG_ALFA
        cout << " a_: " << a_ << " b1: " << b1 << " b2: "<< b2 << " ]]" << endl;
#endif

        th1 = min(RMSE,th1);
      }
      else
      {
#ifdef DEBUG_ALFA      
        cout << "[good fit\tRMSE: " << RMSE << "\tm_RMSE: " << mabs(m_RMSE);
#endif
        // add a part of points to line for cheking in next step
        b1 = a_;
        a_ = round((b1 + b2)/2.0);
        // update thershold	
        th1 = (RMSE + th1)/2.0;
        RMSE_old = RMSE;
	
#ifdef DEBUG_ALFA	
        cout << "\ta_: " << a_ << "\tb1: " << b1 << "\tb2: "<< b2 << " ]" << endl;		
#endif	
      }
    }
    
    cout << "\n" << a_ << "  - " << k_min << "  = " << a_-k_min << "  >  " << th2 << endl;
// check lenght of a line    
    if ((a_ - k_min) > th2)
    {
      ALFA_line_count++;
      cout << "\033[22;36mALFA_line_count:" << ALFA_line_count << " line length: " << a_ - k_min << "\033[22;39m" << endl;

      if ( ALFA_line_count < 8 )
      {
         // Enabale this line if you want see lines in diffrent color      
        glColor3f(seg_color[ALFA_line_count-1].r,seg_color[ALFA_line_count-1].g,seg_color[ALFA_line_count-1].b);
//       glColor3f(seg_color[1].r,seg_color[1].g,seg_color[1].b);
      }
      else
       glColor3f(seg_color[1].r,seg_color[1].g,seg_color[1].b);
       
      // Draw lines that boundary search determine
      glLineWidth(1);
      glBegin(GL_LINES);
        glVertex2f(seg.p[k_min].x/1000.0-4.0,seg.p[k_min].y/1000.0-4.0);
        glVertex2f(seg.p[a_-1].x/1000.0-4.0,seg.p[a_-1].y/1000.0-4.0);
      glEnd();
      glLineWidth(1);      
      
      // Draw two points in begin and end of the line
      glColor3f(1.0,0.27,0.0);
      glPointSize(3);
      glBegin(GL_POINTS);
        glVertex2f(seg.p[k_min].x/1000.0-4.0,seg.p[k_min].y/1000.0-4.0);
        glVertex2f(seg.p[a_-1].x/1000.0-4.0,seg.p[a_-1].y/1000.0-4.0);
      glEnd();      
      glPointSize(1);
      
      glutSwapBuffers();
      
    }
    k_min = a_ + 2;
  }
}

void ALFA(vector<Point> &p)
{
  for (int j=0; j < segment_boundary.size() ; j++ )
  {
    cout << "\n\n\033[22;32m ---------- Segment " << j+1 <<" ---------\n\n\033[22;39m"<< endl;  
    
    int S_len = segment_boundary[j].end-segment_boundary[j].begin+1;
    cout << "Seg Lenght: " << segment_boundary[j].end << " - " << segment_boundary[j].begin << " + 1 = " << S_len << " == " << segment_boundary[j].p.size() << endl << endl;
    
    if ( segment_boundary[j].p.size() < th2)
    {
      cout << segment_boundary[j].p.size() << " Small Segment --> continue." << endl;
      continue;
    }
    
#ifdef DEBUG_SHOW_S
     for (int i=0 ; i <= segment_boundary[j].p.size()-1 ; i++)
     {
       cout << "Index == " << i << "       X == " << segment_boundary[j].p[i].x << "\tY == " << segment_boundary[j].p[i].y << endl;
     }
#endif    
    
    Sx[0]=0.0;
    Sy[0]=0.0;
    Sxy[0]=0.0;
    Sxx[0]=0.0;
    Syy[0]=0.0;
    for (int i=1; i < segment_boundary[j].p.size() ; i++)
    {
      Sx[i] = Sx[i-1] + segment_boundary[j].p[i-1].x;
      Sy[i] = Sy[i-1] + segment_boundary[j].p[+i-1].y;
      
      Sxy[i] = Sxy[i-1] + segment_boundary[j].p[i-1].x * segment_boundary[j].p[i-1].y;
      Sxx[i] = Sxx[i-1] + segment_boundary[j].p[i-1].x * segment_boundary[j].p[i-1].x;
      Syy[i] = Syy[i-1] + segment_boundary[j].p[i-1].y * segment_boundary[j].p[i-1].y;
#ifdef DEBUG_SHOW_S
      cout << i << " == "<< Sx[i] << "  "  << Sy[i] << "  " << Sxy[i] << "  " << Sxx[i] << "  " <<  Syy[i] << endl;
#endif
    }
    BoundarySearch(segment_boundary[j]);
  }
  ALFA_line_count=0;
}

void SegmentationFinal(const vector<Point> &pp)
{
  cout << "\n\n\033[22;35m-------------------- New Segment - Scan Count: " << scan_count<< " --------------------\033[22;31m\n\n";
  seg_bound segment;
  segment.begin=0;
  vector<Point> vecPTmp;
  Point pTmp;
  for (int i = 0 ; i < pp.size(); i++)
  {
      if (pp[i].r > 0.0)
      {
          pTmp.r = pp[i].r;
          pTmp.x = pp[i].x;
          pTmp.y = pp[i].y;
          pTmp.index = i;
          vecPTmp.push_back(pTmp);
      }
  }
  for (int i=1; i <= vecPTmp.size();i++)
  {
    segment.p.push_back(vecPTmp[i-1]);
    if (mabs(vecPTmp[i].r-vecPTmp[i-1].r) > 200 /*polar step th*/|| mabs(vecPTmp[i].index-vecPTmp[i-1].index) > 5/*ray step th*/ || i == vecPTmp.size() )
    {
      segment.end = i-1;
      segment_boundary.push_back(segment);
      segment.p.clear();
      cout << segment_boundary.size() << ": " << segment.begin << "  " << segment.end << " = " << segment.end-segment.begin+1 << endl;
      segment.begin=i;
    }
  }

  vector<Point> out;
  out.resize(pp.size());
  for (int i = 0 ; i < segment_boundary.size() ; i++)
      if (segment_boundary[i].end - segment_boundary[i].begin > 20)
          for (int j=0; j < segment_boundary[i].p.size() ;j++)
              out[segment_boundary[i].p[j].index] = segment_boundary[i].p[j];
}

void Segmentation1(const vector<Point> &pp)
{
  float offset = 0.0;
  cout << "\n\n\033[22;35m-------------------- New Segment - Scan Count: " << scan_count<< " --------------------\033[22;31m\n\n";
  seg_bound segment;
  segment.begin=0;
  vector<Point> vecPTmp;
  Point pTmp;
  for (int i = 0 ; i < pp.size(); i++)
  {
      if (pp[i].r > 0.0)
      {
          pTmp.r = pp[i].r;
          pTmp.x = pp[i].x;
          pTmp.y = pp[i].y;
          pTmp.index = i;
          vecPTmp.push_back(pTmp);
      }
  }
  for (int i=1; i <= vecPTmp.size();i++)
  {
    segment.p.push_back(vecPTmp[i-1]);
    if (mabs(vecPTmp[i].r-vecPTmp[i-1].r) > 200 /*polar step th*/|| mabs(vecPTmp[i].index-vecPTmp[i-1].index) > 5/*ray step th*/ || i == vecPTmp.size() )
    {
      segment.end = i-1;
      segment_boundary.push_back(segment);
      segment.p.clear();
      cout << segment_boundary.size() << ": " << segment.begin << "  " << segment.end << " = " << segment.end-segment.begin+1 << endl;
      
      if (segment.end-segment.begin+1 > 20) // seg len th
      {
/*        glColor3f(1.2,0.0,0.0);
        glPointSize(4);
        glBegin(GL_POINTS);
          glVertex2f(vecPTmp[segment.end].x/1000.0-0.0, vecPTmp[segment.end].y/1000.0);
        glEnd();*/
      
        if (polar_view)
        {
            glColor3f(0.0,0.0,1.0);
            offset = i*0.3515625/100.0;
            glBegin(GL_POINTS);
                glVertex2f(offset-0.0,vecPTmp[segment.end].r/1000.0-4.0);
            glEnd();
        }
      }
      else
      {
/*        glColor3f(1.2,1.0,0.0);
        glPointSize(4);
        glBegin(GL_POINTS);
          glVertex2f(vecPTmp[segment.end].x/1000.0-0.0, vecPTmp[segment.end].y/1000.0);
        glEnd();*/
      
        if (polar_view)
        {
            glColor3f(1.2,1.2,0.8);
            offset = i*0.3515625/100.0;
            glBegin(GL_POINTS);
                glVertex2f(offset-0.0,vecPTmp[segment.end].r/1000.0-4.0);
            glEnd();
        }
      }
    
      glPointSize(1);
      segment.begin=i;
    }
  }

  vector<Point> out;
  out.resize(pp.size());

  for (int i = 0 ; i < segment_boundary.size() ; i++)
  {
       if (segment_boundary[i].end - segment_boundary[i].begin > 20)
       {
           cout << "Seg: " << i << endl;
           glColor3f(seg_color[i].r,seg_color[i].g,seg_color[i].b);

           for (int j=0; j < segment_boundary[i].p.size() ;j++)
           {
               glBegin(GL_POINTS);
                   glVertex2f(segment_boundary[i].p[j].x/1000.0-0.0, segment_boundary[i].p[j].y/1000.0-0.0);
               glEnd();
               out[segment_boundary[i].p[j].index] = segment_boundary[i].p[j];
           }
       }
  }

  for (int i = 0 ; i < out.size() ; i++)
  {
//      cout << "[" << out[i].x << " " << out[i].y << "]" << endl;
      glColor3f(1.0,1.0,1.0);
      glBegin(GL_POINTS);
          if (out[i].x != 0 || out[i].y != 0)
              glVertex2f(out[i].x/1000.0-4.0, out[i].y/1000.0-0.0);
      glEnd();
  }
}

void Segmentation(vector<Point> &p)
{
  float offset = 0.0;
  cout << "\n\n\033[22;35m-------------------- New Segment - Scan Count: " << scan_count<< " --------------------\033[22;31m\n\n";
  seg_bound segment;
  segment.begin=0;
  for (int i=1; i <= p.size();i++)
  {
    segment.p.push_back(p[i-1]);
    if (mabs(p[i].r-p[i-1].r) > 200 /*polar step th*/|| mabs(p[i].index-p[i-1].index) > 5/*ray step th*/ || i == p.size() )
    {
      segment.end = i-1;
      segment_boundary.push_back(segment);
      segment.p.clear();
      cout << segment_boundary.size() << "- " << segment.begin << "  " << segment.end << " = " << segment.end-segment.begin+1 << endl;
      
      if (segment.end-segment.begin+1 > 20) // seg len th
      {
        glColor3f(1.2,0.0,0.0);
        glPointSize(4);
        glBegin(GL_POINTS);
          glVertex2f(p[segment.end].x/1000.0-0.0, p[segment.end].y/1000.0);
        glEnd();
      
        if (polar_view)
        {
            glColor3f(0.0,0.0,1.0);
            offset = i*0.3515625/100.0;
            glBegin(GL_POINTS);
                glVertex2f(offset-0.0,p[segment.end].r/1000.0-4.0);
            glEnd();
        }
      }
      else
      {
        glColor3f(1.2,1.0,0.0);
        glPointSize(4);
        glBegin(GL_POINTS);
          glVertex2f(p[segment.end].x/1000.0-0.0, p[segment.end].y/1000.0);
        glEnd();
      
        if (polar_view)
        {
            glColor3f(1.2,1.2,0.8);
            offset = i*0.3515625/100.0;
            glBegin(GL_POINTS);
                glVertex2f(offset-0.0,p[segment.end].r/1000.0-4.0);
            glEnd();
        }
      }
    
      glPointSize(1);
      segment.begin=i;
    }
  }
}

vector<Point> get_laser_data_new(float *x, float *y, float *t,bool PolarSegDraw)
{
     float r;
     vector<Point> p;
     float offset = 0.0;
     Point temp;

     if ( feof(f) )
     {
        printf("End Of File\n");
        exit(0);
        return p;
     }
     
     /* Reads odometry data: x,y,and angle, but they aren't using in this program. */
     fscanf(f,"%f %f %f",x,y,t);
#ifdef DEBUG
     printf("------= %f\t%f\t%f\n",*x,*y,*t);
#endif

     if (PolarSegDraw)
       glColor3f(0.0,1.0,0.5);

     /* Reads a entire line of the file. */
     for (int i=0;i<(769-43-44);i++)
     {
       fscanf(f,"%f",&r);
       r /= 1;
       
       // Presents a rnage in polar view without any limitation. Thus it's useful for error monitoring of the ranges.
       if (PolarSegDraw /*&& r > 20.0 && r < 4400.0*/)
       {
          if (/*polar_view*/ r > 0.0)
          {
              glBegin(GL_POINTS);
                  glVertex2f(offset-0.0,r/1000.0-6.0);
              glEnd();
              offset += 0.3515625/100.0;
          }
       }

       /* Sslects usable ranges. */
       if ( true /*r > 0.0 && r < 4400.0*/ )
       {
          temp.x=r*cos((0.3515625*i-30.5859375)*PI/180); //-135+16.1719   // 16.1719=44*0.3515625
          temp.y=r*sin((0.3515625*i-30.5859375)*PI/180);
          temp.r=r;
          temp.index=i;
#ifdef DEBUG
         printf("%d = ray: %f\t%f\t%f\n",i,temp.x,temp.y,r);
#endif
         p.push_back(temp);
       }
     }

     // Draws maximum boundary of laser range finder for Hokuyo. (4.4m)
     if ( PolarSegDraw )
     {
        // 30.5859375 = (360/1024)*(44+43)
        // 30.5859375 = (270/768)*(44+43)
        Point temp1;
        temp.x=4.4*cos((0.3515625*0-30.5859375)*PI/180); //-135+16.1719   // 16.1719=44*0.3515625
        temp.y=4.4*sin((0.3515625*0-30.5859375)*PI/180);

        temp1.x=4.4*cos((0.3515625*682-30.5859375)*PI/180); //-135+16.1719   // 16.1719=44*0.3515625
        temp1.y=4.4*sin((0.3515625*682-30.5859375)*PI/180);

        glColor3f(0.7,0.2,0.3);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        gluDisk(diskQuadric,4.4,4.4,30,1);
        glBegin(GL_LINES);
          glVertex2f(0-0.0,0.0-0.0); glVertex2f(temp.x-0.0,temp.y-0.0);
          glVertex2f(0-0.0,0.0-0.0); glVertex2f(temp1.x-0.0,temp1.y-0.0);
        glEnd();
        glutSwapBuffers();
     }
     return p;
}

vector<Point> get_laser_data_newVersion(float *x, float *y, float *t)
{
     float r,r_temp;
     fscanf(f, "%f %f %f", x , y, t);
#ifdef DEBUG            
     printf("Position: %f %f %f\n", *x, *y, *t);
#endif

     for (int i=0 ; i<0 ; i++)
     {
         fscanf(f,"%f",&r_temp);
#ifdef DEBUG
         printf("%f ",r_temp);
#endif
     }	
#ifdef DEBUG
     printf("\n");
#endif
     
     vector<Point> p1;
     Point temp;
     for ( int i=0; i < (769-43-44);i++)
     {
       fscanf(f,"%f",&r);
/*       char chTMP;
       fscanf(f,"%c",&chTMP);
       if (chTMP == '\n')
       {
           cout << i << " data point. "<< "end of line." << endl;
           break;
       }*/
//       if (r<20)
//        continue;
       if (r == 0)
           continue;
       r *=1000;
       temp.x=r*cos((0.3515625*i-120)*PI/180); //-135+16.1719   // 16.1719=44*0.3515625   
       temp.y=r*sin((0.3515625*i-120)*PI/180);
#ifdef DEBUG       
       printf("%f	%f	%f\n",temp.x,temp.y,r);
#endif
       p1.push_back(temp);
     }     
     
     for (int i=0;i<0;i++)
     {
 	fscanf(f,"%f",&r_temp);
#ifdef DEBUG
 	printf("%f ",r_temp);
#endif
     }	
#ifdef DEBUG
     printf("\n");
#endif     
     
     return p1;
}

vector<Point> get_nao_data(bool PolarSegDraw)
{
     vector<Point> p;
     Point temp;
     float x,y;
     char EOL;

     if ( feof(f) )
     {
        printf("End Of File\n");
        exit(0);
        return p;
     }
     
     /* Reads a entire line of the file. */
     for (int i=0; ;i++)
     {
          fscanf(f,"%f %f%c",&x,&y,&EOL);

          if (x == -100 && y == -100)
          {
              break;
          }

          temp.x = x;
          temp.y = y;
          temp.r = 0;
          temp.index = i;
          temp.hasCluster = false;
          temp.clusterIndex = 0;
#ifdef DEBUG
         printf("%f\t%f\t%d\n",temp.x,temp.y,temp.index);
#endif
         p.push_back(temp);
         
         if (PolarSegDraw)
         {
            glColor3f(0.1,0.2,0.1);
    	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//            gluDisk(diskQuadric,4.4,4.4,30,1);
            
            glColor3f(0.0,0.0,1.0);
            glPointSize(1);
    	    glBegin(GL_POINTS);
                glVertex2f(temp.x-0.0,temp.y-0.0);
    	    glEnd();
//    	    glutSwapBuffers();
//            usleep(20);
         }         

//        if ( EOL == '\n' )
//            break;

        if ( feof(f) )
        {
          printf("End Of File\n");
          exit(0);
          return p;
        }
     }
     return p;
}

//static bool xorder (const Point& p, const Point& q)
//{
//  return p.x < q.x;
//}

static bool yorder (const Point& p, const Point& q) 
{
  return p.y < q.y;
}

int dist( Point a , Point b)
{
  return sqrt(pow(b.x-a.x,2)+pow(b.y-a.y,2));
}

//------ ANN --------
/*bool readPt(Point &a, ANNpoint p)	// read point
{
    p[0] = a.x;
    p[1] = a.y;
    return true;
}*/

//------ ANN --------
/*void printPt(ostream &out, ANNpoint p)	// print point
{
        out << "(" << p[0];
        for (int i = 1; i < dim; i++) {
                out << ", " << p[i];
        }
        out << ")\n";
}*/

/*int RBNN ( vector<Point> &p, float r, int nMin)
{
    //------ ANN --------
    int                 nPts;       // actual number of data points
    ANNpointArray	dataPts;    // data points
    ANNpoint		queryPt;    // query point
    ANNidxArray		nnIdx;      // near neighbor indices
    ANNdistArray	dists;      // near neighbor distances
    ANNkd_tree* 	kdTree;     // search structure

    queryPt = annAllocPt(dim);           // allocate query point
    dataPts = annAllocPts(maxPts, dim);  // allocate data points

    nnIdx = new ANNidx[k];		// allocate near neigh indices
    dists = new ANNdist[k];		// allocate near neighbor dists

    nPts = 0;				// read data points
    //------ ANN --------


    int numofCluster = 0;

    while ( nPts < maxPts && nPts < p.size() )
    {
        readPt(p[nPts], dataPts[nPts]);
#ifdef DEBUG_RBNN
        printPt(cout, dataPts[nPts]);
#endif
        nPts++;
    }

    kdTree = new ANNkd_tree(				// build search structure
                                    dataPts,		// the data points
                                    nPts,		// number of points
                                    dim);		// dimension of space
    for (int i=0 ; i < p.size() ; i++)
    {
//        if ( p[i].hasClprintuster )
//            continue;

        queryPt = dataPts[i];

        kdTree->annkSearch(			// search
                        queryPt,		// query point
                        k,			// number of near neighbors
                        nnIdx,			// nearest neighbors (returned)
                        dists,			// distance (returned)
                        eps);			// error bounds

#ifdef DEBUG_RBNN
        cout << endl << "dists: ";
        for (int m = 0 ; m < k ; m++)
            cout << dists[m] << '\t';
        cout << endl;
#endif

        for ( int j = 0; j < k ; j++)
        {
            if ( p[i].hasCluster && p[nnIdx[j]].hasCluster &&
                 p[i].clusterIndex == p[nnIdx[j]].clusterIndex )
            {
#ifdef DEBUG_RBNN
                cout << "-- Contunue --" << endl;
#endif
                continue;
            }
            else if ( p[i].hasCluster && p[nnIdx[j]].hasCluster && dists[j] < 0.0005 )
            {
#ifdef DEBUG_RBNN
                cout << "-- Merge --" << p[i].clusterIndex << " , " << p[nnIdx[j]].clusterIndex << "\tDist: " << dists[j] << endl;
#endif
            }
            else if ( (p[i].hasCluster || p[nnIdx[j]].hasCluster) && dists[j] < 0.5 )
            {
                if (p[i].hasCluster)
                {
                    p[nnIdx[j]].hasCluster = true;
                    p[nnIdx[j]].clusterIndex = p[i].clusterIndex;
                    segment_boundary[p[i].clusterIndex-1].p.push_back(p[nnIdx[j]]);
                }
                else
                {
                    p[i].hasCluster = true;
                    p[i].clusterIndex = p[nnIdx[j]].clusterIndex;
                    segment_boundary[p[nnIdx[j]].clusterIndex-1].p.push_back(p[i]);
                }
#ifdef DEBUG_RBNN
                cout << "--Add--" << "\tDist: " << dists[j] << endl;
#endif
            }
            else if ( !p[i].hasCluster && !p[nnIdx[j]].hasCluster && dists[j] < 0.005 )
            {
                numofCluster++;
                p[nnIdx[j]].hasCluster = true;
                p[nnIdx[j]].clusterIndex = numofCluster;

                p[i].hasCluster = true;
                p[i].clusterIndex = numofCluster;

                seg_bound temp_seg;

                temp_seg.p.push_back(p[nnIdx[j]]);
                temp_seg.p.push_back(p[i]);

                segment_boundary.push_back(temp_seg);

#ifdef DEBUG_RBNN
                cout << "--New--" << "\tDist: " << dists[j] << endl;
#endif
            }
            else
            {
                cout << "--Out-lier--" << endl;
            }
        }
    }

    cout << "\n\tNN:\tIndex\tDistance\n";
    for (int j=0; j < k ; j++)
    {
       dists[j] = sqrt(dists[j]);		// unsquare distance
       cout << "\t" << j << "\t" << nnIdx[j] << "\t" << dists[j] << "\n";

       glColor3f(seg_color[0].r,seg_color[0].g,seg_color[0].b);
       glBegin(GL_LINES);
//           glVertex2f( queryPt[0], queryPt[1]);
//           glVertex2f( dataPts[nnIdx[j]][0], dataPts[nnIdx[j]][1] );
       glEnd();
    }

    //------ ANN --------
    delete [] nnIdx;	// clean things up
    delete [] dists;
    delete kdTree;
    annClose();		// done with ANN
    //------ ANN --------

#ifdef DEBUG_RBNN
    for (int i=0; i < p.size(); i++)
    {
        cout << "\nCluster Index Point[" << i << "]: " << p[i].clusterIndex << endl;
    }
#endif

    cout << "\nNumber of Cluster: " << numofCluster << endl;
    return numofCluster;
}*/

void defualtDraw( int step )
{
   // Prepares glut scene
   glLoadIdentity();
   glTranslatef(0.0,0.0,transZ);
   glColor3f(0.1,0.1,0.3);

   glBegin(GL_LINE);
     glVertex2f(-20,0);
     glVertex2f(20,0);
     glVertex2f(0,-20);
     glVertex2f(0,20);
   glEnd();
   
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   for (int i=1;i<10;i++)
      gluDisk(diskQuadric,i,i,30,3);
   //------ glut ------

   vector<Point> p;

   /* Reads a line of the file */
//   p = get_nao_data(true);
   float xx1,yy1,tt1,xx2,yy2,tt2;
   p=get_laser_data_new(&xx1,&yy1,&tt1,true);

   Segmentation1(p);

/*   if ( p.size() > k )
       RBNN(p,10,10);

  for (int i=0; i < p.size() ;i++)
  {
      if ( p[i].clusterIndex <= 11 )
          glColor3f(seg_color[p[i].clusterIndex-1].r,seg_color[p[i].clusterIndex-1].g,seg_color[p[i].clusterIndex-1].b);
      else
          glColor3f(1.0, 0.5, 0.7); // Use same color for other Cluster
      glPointSize(3);
      glBegin(GL_POINTS);
          glVertex2f(p[i].x-5.0, p[i].y);
      glEnd();
    glutSwapBuffers();
  }*/

   /* Finds out the line boundaries. */
   ALFA(p);

// Enable this section to draw the segmentation with various color.
/*   for (int j=0; j < segment_boundary.size() ; j++)
   {
     glColor3f(seg_color[j].r,seg_color[j].g,seg_color[j].b);

     for (int i=0; i < segment_boundary[j].p.size() ;i++)
     {
       glBegin(GL_POINTS);
         glVertex2f(segment_boundary[j].p[i].x/1000.0-0.0, segment_boundary[j].p[i].y/1000.0-0.0);
       glEnd();
     }
   }*/

   segment_boundary.clear();
//   printf("Shot Index = %d\n",++scan_count);

   /* Draws current line of the log file (Laser Range Finder data) on the cartesian plane. */
/*   glColor3f(1.0, 0.0, 1.0);
   for (int i=0; i < p.size();i++)
   {
      glBegin(GL_POINTS);
        glVertex2f(p[i].x/1000.0+0.0, p[i].y/1000.0+0.0);
      glEnd();
   }*/
   glutSwapBuffers();
}

void OnKey(unsigned char key, int x , int y)
{
   switch ( key )
   {
      /* It clears the scene and skips to the next 10th step on the log file. */
      case 'm':
      {
         glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);   
         defualtDraw( 10 );
      }break;
      /* It clears the scene and moves to the next step on the log file. */
      case 'c':
      {
         glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
         defualtDraw( 1 );
      }break;
      /* It clears just the scene. */
      case 'd':
      {
         glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);   
      }break;
      /* It clears and zooms in the scene, and moves to the next step on the log file. */
      case 'a':
      {
         transZ -= 1.0;
         glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);   
         defualtDraw( 1 );
      }break;
      /* It clears and zooms out the scene, and moves to the next step on the log file. */
      case 'z':
      {
         transZ += 1.0;
         glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
         defualtDraw( 1 );
      }break;
      /* Reserved, it act same as case 'c' now. */
      case 'l' :
      {
         glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);   
         defualtDraw( 1 );
      }break;
      /* Quit the program. */
      case 'q' :
      {
        exit(1);
      }break;
   }
}

void OnReshape(int w, int h)
{
   if ( h==0 )
      h=1;
   
   // set the drawable region of the window
   glViewport(0,0,w,h);
   
   // set up the projection matrix 
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   
   // just use a perspective projection
   gluPerspective(45.0,w/h,1.0,500);
   
   // go back to modelview matrix so we can move the objects about
   glMatrixMode(GL_MODELVIEW);
   
   glClearColor(0.0,0.0,0.0,0.0);
   
   // clear the screen & depth buffer
   glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
   
   defualtDraw( 0 );
}

void OnDraw() 
{
   glutSwapBuffers();
}

void OnInit() 
{
   diskQuadric = gluNewQuadric();
}

void OnExit() 
{
    cout << endl << ":-) Enjoy!" << endl;
}

void OnIdle() {
   // redraw the screen
   glutPostRedisplay();
   
   // Enable this section for automatic processing a log file
//   glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
//   defualtDraw( 1 );
//   usleep(90000);
}

int main( int argc, char* argv[] )
{
   seg_color[0].r=  0.0; seg_color[0].g= 0.0;  seg_color[0].b= 1.0; // Blue
   seg_color[1].r=  0.0; seg_color[1].g= 1.0;  seg_color[1].b= 0.0; // Green
   seg_color[2].r=  0.0; seg_color[2].g= 1.0;  seg_color[2].b= 1.0; // Cyan
   seg_color[3].r=  1.0; seg_color[3].g= 0.0;  seg_color[3].b= 0.0; // Red
   seg_color[4].r=  1.0; seg_color[4].g= 0.0;  seg_color[4].b= 1.0;
   seg_color[5].r=  1.0; seg_color[5].g= 1.0;  seg_color[5].b= 0.0;
   seg_color[6].r=  1.0; seg_color[6].g= 1.0;  seg_color[6].b= 1.0; // White
   seg_color[7].r=  1.0; seg_color[7].g= 0.27; seg_color[7].b= 0.0;
   seg_color[8].r=  1.0; seg_color[8].g= 0.07; seg_color[8].b= 0.57;
   seg_color[9].r=  0.0; seg_color[9].g= 0.5;  seg_color[9].b= 0.6; // Beutiful Blue
   seg_color[10].r= 0.0; seg_color[10].g=0.3;  seg_color[10].b=0.0; // Dark Green
   
   if (argc < 3)
   {
      printf("\nPlease enter address of log file.\n for instance: ./line ./nao.log Scene_trans\n\n");
      exit(0);
   }
   else if (argc == 3 )
   {
     f=fopen(argv[1],"r");
     transZ = atoi(argv[2]);

     /* Skip lines of the file*/
//     vector<Point> temp;
//     for (int i=0; i< 2 ; i++)
//       temp = get_nao_data(false);
   }
   
   // initialise glut
   glutInit(&argc,argv);
   
   // request a depth buffer, RGBA display mode, and we want double buffering
   glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);
   
   // set the initial window size
   glutInitWindowSize(800,800);
   
   // create the window
   glutCreateWindow("MRL-SPL Landmark Extraction Power by ALE");
   
   // set the function to use to draw our scene
   glutDisplayFunc(OnDraw);
   
   // set the function to use to Keyboard
   glutKeyboardFunc(OnKey);
   
   // set the function to handle changes in screen size
   glutReshapeFunc(OnReshape);
   
   // set the idle callback
   glutIdleFunc(OnIdle);
   
   // run our custom initialisation
   OnInit();
   
   // set the function to be called when we exit
   atexit(OnExit);
   
   // this function runs a while loop to keep the program running.
   glutMainLoop();
   return 0;
}
