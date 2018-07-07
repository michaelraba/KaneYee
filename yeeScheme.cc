#include<iostream>
#include<iomanip>
#include<vector>
#include<cassert>
#include<fstream>
#include<sstream>
#include<cmath>

using namespace std;

int gh=1;
int IE = 4+2*gh-1; // really this is grid and boundary data. //ie,je first number was 5, changing to 200...
int JE = 4+2*gh-1; // does't need to be a global var.

// define Pi=3.1415
constexpr double pi() { return std::atan(1)*4; }

class Matrix
{
  friend class YeeScheme;
public:
  Matrix(int sizeX, int sizeY);
  int GetSizeX() const{return dx;}
  int GetSizeY() const{return dy;}
  double &Element(int x, int y);
  void Print() const;
  void PrintToFile(int tstep) ;
  void makeIdentity() ;

  friend ostream &operator<<(ostream &out, const Matrix &m);
  friend Matrix operator+(Matrix m1, Matrix m2);
  friend Matrix operator*(int fac, Matrix m1);
  friend Matrix operator*(Matrix m1, int fac);
  friend Matrix operator*(Matrix m1, Matrix m2);
  double &operator()(int x, int y);

private:
  vector<vector<double> > data;
  int dx, dy;
};

Matrix::Matrix(int sizeX, int sizeY) : dx(sizeX), dy(sizeY)
{ 
  assert(sizeX > 0 && sizeY > 0);
  for(auto i=0; i<dx; ++i)
    {
      vector<double> temp;
      // add for i loop here..
      for(auto j =0; j < dy; ++j)
        temp.push_back(0);
      data.push_back(temp);
    }
}

double &Matrix::operator()(int x,int y)
{
  assert(x >= 0 && x < dx && y >= 0 && y < dy);
  return this->data[x][y];
}

void Matrix::makeIdentity() 
{
    for (  auto x =0; x < dx;++x)
      for (auto y =0; y < dy;++y){
      this->data[x][y] = 1;
  }
}

//model function that works with correct constructor
void Matrix::Print() const
{
  cout <<endl;
  for(auto x =0; x < dx; ++x){
    for(auto y =0; y < dy; ++y)
      cout << setprecision(3)<< scientific << data[x][y]<< "\t";
    cout<<endl;
  } }


class YeeScheme{
friend class Matrix;
public:
  YeeScheme();
  YeeScheme(vector<Matrix> m);
  void Print();
  /* Obsolete methods for updating the interior
   * These are replaced by the single method YeeScheme::updateInterior()

  vector<Matrix> & updateDz(); //this data structure is a problem because only applu to TE.
  vector<Matrix> & updateEz();
  vector<Matrix> & updateHx();
  vector<Matrix> & updateHy();
  */
  enum PulseOptions { GAUSS_SOURCE, SINE_SOURCE /*...*/ };
  vector<Matrix> & updatePulse(int tStep, PulseOptions);
  enum BoundaryOptions { NEUMANN_BOUNDARY, DIRICHLET_BOUNDARY, TEST, BERENGER_PML /*...*/ };
  enum ModeOptions { TM_MODE, TE_MODE, TEST_CASE /*...*/}; //maybe this enum belongs in main!
  vector<Matrix> & updateInterior(ModeOptions, vector<Matrix>::iterator it , int tStep);
  vector<Matrix> & iterateSolution(int tStep, ModeOptions);
  vector<Matrix> & updateBoundary(BoundaryOptions, vector<Matrix>::iterator it);
  // the first arg of updateboundary should be a vector!
  // each entry of the vector is for the sides.
  vector<BoundaryOptions> FourSidedBoundary;
  vector<Matrix> & updateBoundary(vector<BoundaryOptions>, vector<Matrix>::iterator it);
  vector<Matrix> & updatePML(BoundaryOptions, vector<Matrix>::iterator it);
private:
    vector<Matrix> field;  // need to initialize these in constructur...
};
// ********************************
// Constructor for YeeScheme Class
// ********************************
YeeScheme::YeeScheme()
{ 
 
}
//
YeeScheme::YeeScheme(vector<Matrix> m)
{ 
  field=m; 
}

ostream &operator<<(ostream &out, const Matrix &m)
{
  for(int x=0; x< m.dx; ++x){
      for(int y=0; y< m.dy; ++y)
      //out << m.data[x][y] << "\t";
      //out<< endl;

      out << setprecision(2)<< scientific << m.data[x][y]<< "\t\t";
    out<<endl;
  }
  return out;       
}

vector<Matrix> & YeeScheme::updateBoundary(YeeScheme::BoundaryOptions b, vector<Matrix>::iterator it)
{
  switch (b)
    {
    case NEUMANN_BOUNDARY:
      {

                       /********  GhostPoint! ********* / 
                                  /*ghostPointLoop->vertical*/
                       for (auto x =0; x <   it->dx ;++x){
                         //for (auto y =1; y < it->dy-1;++y){
                         it->data[x][0] = it->data[x][1]; // 88888; //ghostpointOben
                           it->data[x][JE-gh] = it->data[x][JE-gh-1]; // 666666; //ghostpointUnten
                       } /* end boundary->vertical*/ 

                       /* ghostPointLoop->horizontal*/
                       //for (auto x =gh; x <   it->dx-1 ;++x){
                       for (auto y =0; y < it->dy;++y){ //do this one for horiz..
                         it->data[0][y] =it->data[1][y] ; //      444444;
                         it->data[JE-gh][y] =it->data[JE-gh-1][y]  ; //555555;
                       } /* end boundary->horizontal*/ 
  return field;
      }

    case DIRICHLET_BOUNDARY:
      {cout<<"";

                       /********  GhostPoint! ********* / 
                                  /*ghostPointLoop->vertical*/
                       for (auto x =0; x <   it->dx ;++x){
                         //for (auto y =1; y < it->dy-1;++y){
                         it->data[x][0] = (-1)*it->data[x][1]; // 88888; //ghostpointOben0;
                         it->data[x][JE-gh] = (-1)*it->data[x][JE-gh-1];;
                       } /* end boundary->vertical*/ 

                       /* ghostPointLoop->horizontal*/
                       //for (auto x =gh; x <   it->dx-1 ;++x){
                       for (auto y =0; y < it->dy;++y){ //do this one for horiz..
                         it->data[0][y] =(-1)*it->data[1][y] ;;
                         it->data[JE-gh][y] =(-1)*it->data[JE-gh-1][y]  ;
                       } /* end boundary->horizontal*/ 

         
                       //cout<< "DirichletBdyWithParameters:\n"<< *it << endl;
        return field;}

    case TEST:
      {

                       /********  GhostPoint! ********* / 
                                  /*ghostPointLoop->vertical*/
                       for (auto x =0; x <   it->dx ;++x){
                         //for (auto y =1; y < it->dy-1;++y){
                         it->data[x][0] =  88888; //ghostpointOben
                           it->data[x][JE-gh] = 666666; //ghostpointUnten
                       } /* end boundary->vertical*/ 

                       /* ghostPointLoop->horizontal*/
                       //for (auto x =gh; x <   it->dx-1 ;++x){
                       for (auto y =0; y < it->dy;++y){ //do this one for horiz..
                         it->data[0][y] =444444;
                         it->data[JE-gh][y] =555555;
                       } /* end boundary->horizontal*/ 

                          cout<< "Oink!\n";
                          it->Print();
        return field;}
    }
  return field;
}

vector<Matrix> & YeeScheme::updateBoundary(vector<YeeScheme::BoundaryOptions> b, vector<Matrix>::iterator it)
{

  //  switch (b)
  //    {
  //    case NEUMANN_BOUNDARY:
  //      {
  //
  //                       /********  GhostPoint! ********* / 
  //                                  /*ghostPointLoop->vertical*/
  //                       for (auto x =0; x <   it->dx ;++x){
  //                         //for (auto y =1; y < it->dy-1;++y){
  //                         it->data[x][0] = it->data[x][1]; // 88888; //ghostpointOben
  //                           it->data[x][JE-gh] = it->data[x][JE-gh-1]; // 666666; //ghostpointUnten
  //                       } /* end boundary->vertical*/ 
  //
  //                       /* ghostPointLoop->horizontal*/
  //                       //for (auto x =gh; x <   it->dx-1 ;++x){
  //                       for (auto y =0; y < it->dy;++y){ //do this one for horiz..
  //                         it->data[0][y] =it->data[1][y] ; //      444444;
  //                         it->data[JE-gh][y] =it->data[JE-gh-1][y]  ; //555555;
  //                       } /* end boundary->horizontal*/ 
  //  return field;
  //      }
  //
  //    case DIRICHLET_BOUNDARY:
  //      {cout<<"";
  //
  //                       /********  GhostPoint! ********* / 
  //                                  /*ghostPointLoop->vertical*/
  //                       for (auto x =0; x <   it->dx ;++x){
  //                         //for (auto y =1; y < it->dy-1;++y){
  //                         it->data[x][0] = (-1)*it->data[x][1]; // 88888; //ghostpointOben0;
  //                         it->data[x][JE-gh] = (-1)*it->data[x][JE-gh-1];;
  //                       } /* end boundary->vertical*/ 
  //
  //                       /* ghostPointLoop->horizontal*/
  //                       //for (auto x =gh; x <   it->dx-1 ;++x){
  //                       for (auto y =0; y < it->dy;++y){ //do this one for horiz..
  //                         it->data[0][y] =(-1)*it->data[1][y] ;;
  //                         it->data[JE-gh][y] =(-1)*it->data[JE-gh-1][y]  ;
  //                       } /* end boundary->horizontal*/ 
  //
  //         
  //                       //cout<< "DirichletBdyWithParameters:\n"<< *it << endl;
  //        return field;}
  //
  //    case TEST:
  //      {
  //
  //                       /********  GhostPoint! ********* / 
  //                                  /*ghostPointLoop->vertical*/
  //                       for (auto x =0; x <   it->dx ;++x){
  //                         //for (auto y =1; y < it->dy-1;++y){
  //                         it->data[x][0] =  88888; //ghostpointOben
  //                           it->data[x][JE-gh] = 666666; //ghostpointUnten
  //                       } /* end boundary->vertical*/ 
  //
  //                       /* ghostPointLoop->horizontal*/
  //                       //for (auto x =gh; x <   it->dx-1 ;++x){
  //                       for (auto y =0; y < it->dy;++y){ //do this one for horiz..
  //                         it->data[0][y] =444444;
  //                         it->data[JE-gh][y] =555555;
  //                       } /* end boundary->horizontal*/ 
  //
  //                          cout<< "Oink!\n";
  //                          it->Print();
  //        return field;}
  //    }
  return field;
}
























vector<Matrix> & YeeScheme::updatePulse(int tStep, PulseOptions p)
{

  switch (p)
    {
    case GAUSS_SOURCE:
      {
  int ic = IE/2;
  int jc = JE/2;
  double spread = 6.0;
  double T = tStep;
  double t0= 20.0;
  double  pulse = 1*exp(-0.5*(pow((t0-T)/spread,2)));
  cout <<"********* PulseValue is = "<< pulse<<"\n"<<endl;
   field[0].data[ic][jc] = pulse;
  //field[0].data[ic][jc] = 777;
      } //end gauss_source case
    case SINE_SOURCE:
      {
        int ic = IE/2;
        int jc = JE/2;
        float ddx = 0.01;
        float dt = ddx/6e8;
        double spread = 6.0;
        double T = tStep;
        double t0= 20.0;
        double  pulse = 2*pi()*1500*1e6*dt*T;
      } //end sine_source case
    }
  return field;
}

void Matrix::PrintToFile(int tstep) 
{
  string num;
  if (tstep < 10)
    //    num = '0' + to_string(tstep);
  num =  to_string(tstep);
  else
    num = to_string(tstep);
    
  ofstream fieldOutput;
  // string myString = "MyField" + num;
  string myString = "./Field/MyField" + num;
  fieldOutput.open(myString);
  for(auto x =0; x < dx; ++x){
    for(auto y =0; y < dy; ++y)
      fieldOutput << setprecision(2)<< scientific << data[x][y]<< " ";
    fieldOutput<<endl;
  }
  fieldOutput.close();
}



Matrix operator+(Matrix m1, Matrix m2)
{
assert(m1.dx == m2.dx && m1.dy ==m2.dy);
Matrix result = m1;
for(auto i=0; i< m1.dx; ++i)
for(auto j=0; j< m1.dy; ++j)
    result(i,j) = m1(i,j) + m2(i,j);
return result;
}




/* Function: updateInterior()
 * parameters:
 ** ModeOptions p
 ** iterator it
 * Returns: the matrix of vector fields

 Description: Updates the interior (ie. non-boundary and not the PML region)
 modeoptions are:
 *** TE_MODE

 Uses the vector field tuple (Hz, Ex, Ey) which is stored in the vectors of matrices vector<Matrix>. The tuple is updated according to the rules as described in Yee

     Numerical solution of initial boundary value problems involving maxwell's equations in isotropic media

     *** TM_MODE

     Uses the vector field tuple (Ez, Hx,Hy) which is stored in the vectors of matrices vector<Matrix>. The tuple is updated according to the rules as described in Yee. 
*/
vector<Matrix> &  YeeScheme::updateInterior(YeeScheme::ModeOptions m, vector<Matrix>::iterator it, int tStep)
{
         switch (m)
           {
             /*
               |---------+----------------------------|
               | Tm Mode |                            |
               |---------+----------------------------|
               | field   | corresponding vector index |
               |---------+----------------------------|
               | dz      | field[0]                   | 
               | hx      | field[1]                   |
               | hy      | field[2]                   |
               | ez      | field[3]                   |
               | ga      | field[4]                   |
               |---------+----------------------------| */
           case TM_MODE: 
             if (std::distance(it, field.begin()) == 0)
               {
                 cout<<"Dz!\n";
             // * Update dz-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){ //fixed for spurious boundary!
                 field[0].data[x][y] +=  + 0.5*(field[2].data[x][y] - field[2].data[x-1][y] - field[1].data[x][y] + field[1].data[x][y-1] );
               }} // end for loop
             // needs to update boundary for each field iteration!
             return field;
}

             //        if (std::distance(it, field.begin()) == 1)
             //          {
             //        // * Update ez-field
             //        for (auto y =gh; y < field[0].dy-1;++y){
             //          for (auto x =gh; x <   field[0].dx-1 ;++x){
             //            // field[4].data[x][y] +=  0 ; /* *todo* */
             //          }} //end ez-field update
             //        return field; 
             //          } // end iterator if

             else if (std::distance(field.begin(),it) == 1)
               {
                 updatePulse(tStep, GAUSS_SOURCE); //just update before hx..
                 cout<<"dzField SecondPrint:\n"<<field[0]<<"\n\n";
                 cout<<"Hx!\n";
             // * Update hx-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[1].data[x][y] +=  + 0.5*(field[0].data[x][y] - field[0].data[x][y+1] );}} // end hx-field update
             return field;
               }


             else if (std::distance(field.begin(),it) == 2)
               {
                 cout<<"Hy!\n";
             // * Update hy-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[2].data[x][y] += + 0.5*(field[0].data[x+1][y] - field[0].data[x][y] ); 
               }}  // end hy-field update
             return field;
               }

             /*
             |---------+----------------------------|
             | Te Mode |                            |
             |---------+----------------------------|
             | field   | corresponding vector index |
             |---------+----------------------------|
             | hz      | field[0]                   |
             | ex      | field[1]                   |
             | ey      | field[2]                   |
             |---------+----------------------------| */
           case TE_MODE:

             if (std::distance(it, field.begin()) == 0)
               {
             // * Update hz-field
                 cout<<"hzFld\n";
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[0].data[x][y] +=.5*(field[2].data[x][y]-field[2].data[x+1][y])+.5*(field[1].data[x][y+1]-field[1].data[x][y])  ;}} // end hz-field update
             return field;
               }

             else if (std::distance(it, field.begin()) == 1)
               {

                 updatePulse(tStep, GAUSS_SOURCE); //just update before ex...
                 cout<<"exFld\n";
             // * Update ex-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[1].data[x][y] += .5*( field[2].data[x][y] - field[2].data[x][y-1] );}} // end ex-field update
             return field;
               }

             else if (std::distance(it, field.begin()) == 2)
               {
                 cout<<"eyFld\n";
             // * Update ey-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[2].data[x][y] += .5*(field[2].data[x][y] - field[2].data[x-1][y])  ;}} // end ey-field update
             return field;
               }


           case TEST_CASE:
             
             // * Update hz-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[0].data[x][y] =666;}} // end hz-field update

             // * Update ex-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[1].data[x][y] =666;}} // end ex-field update

             // * Update ey-field
               for (auto x =gh; x <   field[0].dx-1 ;++x){
                 for (auto y =gh; y < field[0].dy-1;++y){
                 field[2].data[x][y] =666;}} // end ey-field update

               return field;
           }
         return field;
}


// fix this...the iterator currently is neceesary because updateboundary uses it.
// buuuut updateinterior does not..
// we still need to implent:
// check what iteration it is..then do only that case..
// with updateboundary !
vector<Matrix> &  YeeScheme::iterateSolution(int tStep, YeeScheme::ModeOptions m)
{ 
  vector<Matrix>::iterator it;
  //  vector<Matrix> field = {dz,hx,hy,ez,ga};
  //  for(auto it = field.begin(); it!=field.end()-1; it++)
  // 
  vector<BoundaryOptions> sides = {YeeScheme::BERENGER_PML,YeeScheme::BERENGER_PML,YeeScheme::BERENGER_PML,YeeScheme::BERENGER_PML};
  for(vector<Matrix>::iterator it = field.begin(); it!=field.end()-1; it++)
    {
  //updateBoundary(YeeScheme::NEUMANN_BOUNDARY);
      updateInterior(YeeScheme::TM_MODE, it, tStep );
      updateBoundary(YeeScheme::DIRICHLET_BOUNDARY, it);
      updateBoundary(sides, it);
      it->PrintToFile(tStep);

  //updateBoundary(YeeScheme::TEST, it);
    }
  return field;
}




//int main(int argc, char* argv[])
//{
//  Matrix foo(2,3);
//  foo.Print(); // this prints
//
//  cout<<"Print Bar:\n";
//  Matrix bar(3,2);
//   bar.Print(); //this does not print
//   bar.makeIdentity();
//   bar.Print(); //this does not print
//  
//  return 0;
//}

int main(int argc, char* argv[])
{

  enum ModeOptions { TM_MODE, TE_MODE/*...*/ }; //keep this enim in YeeScheme (for now..)
    int l,n,i,j,ic,jc,nsteps,npml;
    float ddx,dt,T,epsz,pi,epsilon,sigma,eaf;
    float t0,spread,pulse;
    T  =0;
    if (argc == 1 )
      nsteps = 2;
    else
      {
        istringstream ss(argv[1]);
        int x;
        if (!(ss >> x))
          cerr << "Invalid number " << argv[1] << '\n';
        nsteps = x;}
    ic = IE/2;
    jc = JE/2;
    ddx = 0.01; // spatial step
    dt = ddx/6e8; // yee cell size // why...because CFL condition..t scales different from x..
    epsz=8.8e-12;
    pi=3.14159;

    cout<<"^^^^^^^^^^^^^Begin Iterates:\n";
    ModeOptions m = TE_MODE;
    YeeScheme::BoundaryOptions b = YeeScheme::BERENGER_PML;
    
    switch (m)
      {
      case TM_MODE:
        {
          Matrix ga(IE,JE),dz(IE,JE),ez(IE,JE) ;
          Matrix hx(IE,JE),hy(IE,JE);
          ga.makeIdentity();
          int dSize = IE;
          //where is ex, ey fields? haha! And Hz, lol.
          vector<Matrix> field = {dz,hx,hy,ez,ga};
          //  vector<Matrix> Efield = {ex,ey,ez};
          YeeScheme yee(field); // how would i call yee without field?
          // i maybe do not create the yee object in main. Ok.
          for(n=1; n<nsteps; ++n)
            {
                      cout<< "*****Iteration step:"<< n<<"******" << endl;
              T= T+1;
              yee.iterateSolution(T, YeeScheme::TM_MODE);
              //field[0].PrintToFile(T);
            }
        } //end tm-mode case

      case TE_MODE:
        {
          switch(b)
            {
          case YeeScheme::BERENGER_PML:
            {
              cout<<"*******hi******";
              int cellSize=4;
              int IEx=IE+cellSize;
              int JEx=JE+cellSize;
              Matrix hzx(IEx,JEx),hzy(IEx,JEx);
              Matrix experiment(IE,cellSize);
              experiment.makeIdentity();
              cout << "my experiment:\n"<<experiment;
            } //end berenger_pml case
            } //end boundary switch
          //
          //          Matrix ga(IE,JE); //auxillary matrices..
          //          Matrix ex(IE,JE),ey(IE,JE),hz(IE,JE);
          //          ga.makeIdentity();
          //          int dSize = IE;
          //          //where is ex, ey fields? haha! And Hz, lol.
          //          vector<Matrix> field = {ex,ey,hz,ga};
          //          YeeScheme yee(field); // how would i call yee without field?
          //          // i maybe do not create the yee object in main. Ok.
          //          for(n=1; n<nsteps; ++n)
          //            {
          //                      cout<< "*****Iteration step:"<< n<<"******" << endl;
          //              T= T+1;
          //              yee.iterateSolution(T, YeeScheme::TE_MODE);
          //              //field[0].PrintToFile(T);
          //            }
        } //end te-mode case
      }
return 0;
}
//eof
