 //Ref: Sathish Kumar D and Rahla Rabia M P
//Input Files square cells from raster data
// ADI scheme
//Working good final 
//multiple outlet specified as outlet file

#include"iostream"
#include<conio.h>
#include"fstream"
#include<math.h>
#include <stdlib.h>

//declare global functions used in linpack

int idamax ( int n, double dx[], int incx ); 
void dscal ( int n, double sa, double x[], int incx );
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );  
double r8_abs ( double x );
void dswap ( int n, double x[], int incx, double y[], int incy );
int dgefa ( double a[], int lda, int n, int ipvt[] );
void dgedi ( double a[], int lda, int n, int ipvt[], double det[], 
  double work[], int job );



// multiplication function prototype
int multimat(double **,double *,double *,int );



  
using namespace std;




int main()
{
    
    
int nelems,
    i,j,k,
    n,m,
    t,t1,
    nodatavalue,
    dsrow,dscol,
    **outlet,
    noutlet,
    max_size;


double  **cells,
        **hcells,
        **hstcells,
        **hwater,
        **P,
        **Q,
        **Q1,
        **Q2,
        **Qout,
        *Dx,
        *Dy,
        *hst,
        *ht,
        *B1,
        *B,
        h1,
        sqso,
        r,
        nm,
        alpha,
        min_h,
        min_sqso,
        K1,
        K2,
        del_t,
        r1,
        r2,
        tr,
        tn,
        rain,
		del_x,
        del_y,
        sy,
        z1,z2,z3,z4;


        

// declaration for matrix inverse -variables
 
    
     double *a; 
     double det[2];
     int info;
     int *ipvt;
     int job;
     double *work;                



//initialize minimum values and others

min_h = 0.0000000001;
min_sqso =0.00000000000000000000001;




                   //********************************//
                   //                                //      
                   //  READING FLOW PARAMETER FILE   //                
                   //                                //      
                   //********************************//             
    
    fstream file_open1("d:\\Mannings\\005\\newflow_param005.txt",ios::in);
    file_open1>>del_t>>r1>>r2>>tr>>tn>>alpha;
   
    file_open1.close();
    
     fstream file_open5("d:\\Mannings\\005\\Mannings.txt",ios::in);
      file_open5>>nm;

                   //********************************//
                   //                                //      
                   //       READING DEM FILE         //                
                   //                                //      
                   //********************************//             
                                        
    fstream file_open2("d:\\Mannings\\005\\dem005.txt",ios::in);
    file_open2>>m>>n>>del_x>>del_y>>nodatavalue;

    cells = (double **)malloc((n+1)*sizeof(double *));
    hwater= (double **)malloc((n+1)*sizeof(double *));
    hcells = (double **)malloc((n+1)*sizeof(double *));
    hstcells = (double **)malloc((n+1)*sizeof(double *));
    Q1 = (double **)malloc((n+1)*sizeof(double *));
    Q2 = (double **)malloc((n+1)*sizeof(double *));
    Qout = (double **)malloc((n+1)*sizeof(double *));
    outlet = (int **)malloc((n+1)*sizeof(int *));
    
    for(i=1;i<=n;i++)
    {
    cells[i] = (double *)malloc((m+1)*sizeof(double));
    hwater[i] = (double *)malloc((m+1)*sizeof(double));
    hcells[i] = (double *)malloc((m+1)*sizeof(double));
    hstcells[i] = (double *)malloc((m+1)*sizeof(double));
    Q1[i] = (double *)malloc((n+1)*sizeof(double));
    Q2[i] = (double *)malloc((n+1)*sizeof(double));
    Qout[i] = (double *)malloc((n+1)*sizeof(double));
    outlet[i] = (int *)malloc((n+1)*sizeof(int));
    }    
       
 
  
  if(m>=n)  max_size =m;
  else max_size =n;   
  
  
  
    
    for(j=1;j<=m;j++)
    for(i=1;i<=n;i++)
    file_open2>>cells[i][j];
    file_open2.close();
    
   fstream file_open3("d:\\Mannings\\outlet.txt",ios::in); 
   for(j=1;j<=m;j++)
   for(i=1;i<=n;i++)
   file_open3>>outlet[i][j];
   file_open3.close(); 
   
   noutlet=0;
   for(j=1;j<=m;j++)
   for(i=1;i<=n;i++)
   if(outlet[i][j]!=0)
   noutlet+=1;  
    
    
    
    for(j=1;j<=m;j++)
    for(i=1;i<=n;i++)
    hcells[i][j]=hstcells[i][j]=nodatavalue;
    
    for(j=1;j<=m;j++)
    for(i=1;i<=n;i++)   
    hcells[i][j]=cells[i][j];
    

   


    r=1.0;
    
    
     
//intialize with small depth
for(i=1;i<=n;i++)
for(j=1;j<=m;j++)
if(cells[i][j]!=nodatavalue)
{
hwater[i][j]=0;
hcells[i][j]=hcells[i][j]+hwater[i][j];
}    

fstream file_open4("d:\\Mannings\\005\\output_fvm_overland1.txt",ios::out); 

t1 = int(tn/del_t);


//Start Time loop

for(t=0;t<t1;t++)

{
if((t*del_t)<=tr) 
rain = r1;
else
rain = r2;


//Memory allocation P[i][i], Dx, Dy, hst, B with dimension i



Dx = (double *)malloc((n+1) * sizeof(double ));
Dy =(double *)malloc((n+1) * sizeof(double ));
hst =(double *)malloc((n+1) * sizeof(double ));
B = (double *)malloc((n+1) * sizeof(double ));
P = (double **)malloc((n+1) * sizeof(double *));
for(i=1;i<=n;i++)
P[i]=(double *)malloc((n+1) * sizeof(double));



for(j=1;j<=m;j++)
{
  
 for(i=1;i<=n;i++)
    {
       Dx[i]=0; 
       Dy[i]=0; 
       hst[i]=0;
       B[i]=0;
       for(k=1;k<=n;k++)
       P[i][k]=0;
    }   
    

 for(i=1;i<=n;i++)
    {
   
  
   if(i!=1)
   {
     if((cells[i][j]!= nodatavalue)&&(cells[i-1][j]==nodatavalue))
     K1=0;    
     else if((cells[i][j]!= nodatavalue)&&(cells[i-1][j]!=nodatavalue))
     {
        h1 = 0.5*(hwater[i-1][j]+hwater[i][j]);
        
		if(j!=m)
			{
			if(cells[i][j+1] == nodatavalue) 
					z1 = hcells[i][j];
			else
					z1 = hcells[i][j+1];
			if(cells[i-1][j+1] == nodatavalue) 
					z2 = hcells[i-1][j];
			else
					z2 = hcells[i-1][j+1];
			}
		else
			{
					z1 = hcells[i][j];
					z2 = hcells[i-1][j];
			}
		if(j!=1)
			{
			if(cells[i][j-1] == nodatavalue) 
					z3 = hcells[i][j];
			else
					z3 = hcells[i][j-1];
			if(cells[i-1][j-1] == nodatavalue) 
					z4 = hcells[i-1][j];
			else
					z4 = hcells[i-1][j-1];
			}
        
		else
			{
					z3 = hcells[i][j];
					z4 = hcells[i-1][j];
			}

        sqso = (1/(del_x*del_x))*(((hcells[i][j]-hcells[i-1][j])*(hcells[i][j]-hcells[i-1][j]))+
               (r*r/16)*((z1+ z2- z3- z4)*(z1+z2-z3-z4)));
     
      if(h1>0)
	  K1= (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
	  else
	  K1=0;


	  if(hcells[i][j]>hcells[i-1][j])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
	  else if(hcells[i-1][j]>hcells[i][j])
		  {
		  if((hwater[i-1][j]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
      else if(hcells[i-1][j]==hcells[i][j])
		  K1=0;
      
      }
    } 
    else
    K1=0;   
     
  
      
     if(i!=n) 
     { 
     if((cells[i][j]!= nodatavalue)&&(cells[i+1][j]==nodatavalue))
      K2=0;
     else if((cells[i][j]!= nodatavalue)&&(cells[i+1][j]!=nodatavalue))
      {
          h1 = 0.5*(hwater[i][j]+hwater[i+1][j]);

		if(j!=m)
		{
          if(cells[i][j+1] == nodatavalue) 
                    z1 = hcells[i][j];
          else      
                    z1 = hcells[i][j+1]; 
		  
          if(cells[i+1][j+1] == nodatavalue) 
                    z2 = hcells[i+1][j];
          else
                    z2 = hcells[i+1][j+1];
		}
		else
		{
				z1 = hcells[i][j];
				z2 = hcells[i+1][j];	
		}

        if(j!=1)
		{
		
		if(cells[i][j-1] == nodatavalue) 
                    z3 = hcells[i][j];
          else
                    z3 = hcells[i][j-1];

          if(cells[i+1][j-1] == nodatavalue) 
                    z4 = hcells[i+1][j];
          else
                    z4 = hcells[i+1][j-1];
		}
          
        else
		{
			z3 = hcells[i][j];
			z4 = hcells[i+1][j];
		}
		
		sqso = (1/(del_x*del_x))*(((hcells[i+1][j]-hcells[i][j])*(hcells[i+1][j]-hcells[i][j]))+
                 (r*r/16)*((z2+ z1- z4- z3)* (z2+ z1- z4- z3)));
      
        if(h1>0)
        K2= (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
        else
        K2 =0;
        
        
        
		if(hcells[i][j]>hcells[i+1][j])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
		else if(hcells[i+1][j]>hcells[i][j])
		  {
		  if((hwater[i+1][j]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }		
		else if(hcells[i+1][j]==hcells[i][j])
				K2=0;
      }                 
    } 
    else 
    K2 =0;   
  

   
        
      if((i!=1)&&(i!=n)&&(cells[i][j] != nodatavalue))
     {
     P[i][i-1]= -(0.5*del_t*K1)/(del_x*del_y);
     P[i][i] = 1+((0.5*del_t*(K1+K2))/(del_x*del_y));
     P[i][i+1]= -(0.5*del_t*K2)/(del_x*del_y);
     }    
     else if((i==1)&&(i!=n)&&(cells[i][j] != nodatavalue))
     {
     P[i][i] = 1+((0.5*del_t*(K1+K2))/(del_x*del_y));
     P[i][i+1]= -(0.5*del_t*K2)/(del_x*del_y);
     }  
     else if((i!=1)&&(i==n)&&(cells[i][j] != nodatavalue)) 
     {
     P[i][i-1]= -(0.5*del_t*K1)/(del_x*del_y);
     P[i][i] = 1+((0.5*del_t*(K1+K2))/(del_x*del_y));
     }
     
     
       
     
         
      if(cells[i][j] != nodatavalue)
      {
       if(K2==0)
       z1=0;
       else 
       z1=K2*(hcells[i+1][j]-hcells[i][j]);
       if(K1==0)
       z2=0;
       else
       z2=K1*(hcells[i-1][j]-hcells[i][j]);
         Dx[i] = ((0.5*del_t)/(del_x*del_y))*(z1+z2);
      }
  
  
  
  if (j!=1)
      {
     if((cells[i][j]!= nodatavalue)&&(cells[i][j-1]==nodatavalue))
     K1=0;
     else  if((cells[i][j]!= nodatavalue)&&(cells[i][j-1]!=nodatavalue))
      {
        h1 = 0.5*(hwater[i][j-1]+hwater[i][j]);
        
		if(i!=n)
			{
			if(cells[i+1][j] == nodatavalue)
					z1= hcells[i][j];
			else 
					z1=hcells[i+1][j];
			if(cells[i+1][j-1] == nodatavalue) 
					z2= hcells[i][j-1];
			else 
					z2=hcells[i+1][j-1];
			}
		else
			{
				z1= hcells[i][j];
				z2= hcells[i][j-1];
			}
        
		if(i!=1)
			{
			if(cells[i-1][j] == nodatavalue) 
					z3= hcells[i][j];
			else
					z3=hcells[i-1][j];
			if(cells[i-1][j-1] == nodatavalue) 
					z4 = hcells[i][j-1];
			else 
					z4= hcells[i-1][j-1];
			}

		else
			{
				z3= hcells[i][j];
				z4 = hcells[i][j-1];
			}
     
        sqso = (1/(del_y*del_y))*(((hcells[i][j]-hcells[i][j-1])*(hcells[i][j]-hcells[i][j-1]))+
               (r*r/16)*((z1+ z2- z3- z4)* (z1+ z2- z3- z4)));
      
      if(h1>0)
      K1= (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
      else
      K1 =0;
      
      if(hcells[i][j]>hcells[i][j-1])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
	  else if(hcells[i][j-1]>hcells[i][j])
		  {
		  if((hwater[i][j-1]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
      else if(hcells[i][j-1]==hcells[i][j])
		  K1=0;

      }
     }    
   else
   K1=0; 
  

 
  if(j!=m)
  {    
      if((cells[i][j]!= nodatavalue)&&(cells[i][j+1]==nodatavalue))
        K2=0;
      else if((cells[i][j]!= nodatavalue)&&(cells[i][j+1]!=nodatavalue))
      {
        h1 = 0.5*(hwater[i][j+1]+hwater[i][j]);
        
		if(i!=n)
			{
			if(cells[i+1][j] == nodatavalue) 
					z1 = hcells[i][j];
			else 
					z1= hcells[i+1][j]; 
			if(cells[i+1][j+1] == nodatavalue)  
					z2= hcells[i][j+1];
			else
					z2=hcells[i+1][j+1];
			}
		else
			{
				z1 = hcells[i][j];
				z2= hcells[i][j+1];
			}
		
		if(i!=1)
			{
			if(cells[i-1][j] == nodatavalue)
					z3 = hcells[i][j];
			else
					z3= hcells[i-1][j];
			if(cells[i-1][j+1] == nodatavalue) 
					z4= hcells[i][j+1];
			else
					z4 =hcells[i-1][j+1];
			}

		else
			{
				z3 = hcells[i][j];
				z4= hcells[i][j+1];
			}
        
        sqso = (1/(del_y*del_y))*(((hcells[i][j+1]-hcells[i][j])*(hcells[i][j+1]-hcells[i][j]))+
               (r*r/16)*((z2+ z1- z4- z3)*(z2+ z1- z4- z3)));
      
      if(h1>0)
      K2 = (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
	  else
	  K2 = 0;
	  if(hcells[i][j]>hcells[i][j+1])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
	  else if(hcells[i][j+1]>hcells[i][j])
		  {
		  if((hwater[i][j+1]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
      else if(hcells[i][j+1]==hcells[i][j])
				K2=0;


      }
  }  
  else
  K2=0;    
     

  

  
      if(cells[i][j] != nodatavalue)
      {
         
         
         if(K2==0)
         z1=0;
         else
         z1= K2*(hcells[i][j+1]-hcells[i][j]);
             
         if(K1==0)
         z2=0;
         else 
         z2 = K1*(hcells[i][j-1]-hcells[i][j]);
         Dy[i] = ((0.5*del_t)/(del_x*del_y))*(z1+z2);
                    
         
   
      }
    
  } //end of i loop
 
  
 
    
     
  for(i=1;i<=n;i++) 
  {
      hst[i]=0;
      if(cells[i][j] != nodatavalue)
      {
     if(outlet[i][j]==1)
     
         {
         h1 = hwater[i][j];
         if(h1>min_h)
         K2 = pow(h1,(5/3.0))/nm;
         else K2=0;
         sy = (cells[i][j-1]- cells[i][j])/del_y;
         Q1[i][j] = (K2*pow(sy,0.5)*del_x);
        
         B[i]=Dy[i]+hcells[i][j]- (K2*pow(sy,0.5)*del_t/(del_y))+ 0.5*rain*del_t;
         }  
      else  
	  B[i]=Dy[i]+hcells[i][j]+0.5*rain*del_t; 
      }     
      else
      B[i]= 0;
  }
//Remove zero in diagonals of P Matrix

//count non zero elements in diagonal
nelems = 0;
for (i=1;i<=n;i++)
if(P[i][i]!=0)
{
nelems = nelems+1;


}
 




//Allocate memory for Q, ht matrix
 
 

if(nelems!=0)

 {

 Q=(double **)malloc((nelems+1) * sizeof(double *));
 ht = (double *)malloc((nelems+1) * sizeof(double ));
 B1 = (double *)malloc((nelems+1) * sizeof(double ));
 a=(double *)malloc((nelems*nelems) * sizeof(double ));
 ipvt = (int *)malloc((nelems) * sizeof(int));
 work = (double *)malloc((nelems) * sizeof(double));
 for (i=1;i<=nelems;i++)
 Q[i]=(double *)malloc((nelems+1) * sizeof(double ));




//initialize  Q matrix
for (i=1;i<=nelems;i++)
{
ht[i]=0;
for (k=1;k<=nelems;k++)
Q[i][k]=0;
}



//Write non zero elements in Q matrix
for (i=1,k=1;i<=n;i++)
{
    if((P[i][i]!=0)&&(i!=1)&&(i!=n))
    {
      Q[k][k]=P[i][i];  
      if(P[i][i-1]!=0) Q[k][k-1]=P[i][i-1]; 
      if(P[i][i+1]!=0) Q[k][k+1]=P[i][i+1];	 		
      B1[k]=B[i];
      k+=1;
    }
    else if((P[i][i]!=0)&&(i==1)&&(i!=n))
    {
      Q[k][k]=P[i][i];  
      if(P[i][i+1]!=0)Q[k][k+1]=P[i][i+1];
      B1[k]=B[i];
      k+=1;
    }    
    else if((P[i][i]!=0)&&(i!=1)&&(i==n))
    {
      Q[k][k]=P[i][i];  
      if(P[i][i-1]!=0)Q[k][k-1]=P[i][i-1];
      B1[k]=B[i];
      k+=1;
    }    
}    






 //Inverse Q matrix

    for(i=1;i<=nelems;i++)
   for(k=1;k<=nelems;k++)
    a[((i-1))*nelems +(k-1)]= Q[i][k];      
     
    info = dgefa ( a, nelems, nelems, ipvt );

    if ( info != 0 )
    {
       cout << "  Error!  The matrix is nearly singular!\n";
       return 0;
    }        
    job = 11;
    dgedi ( a, nelems, nelems, ipvt, det, work, job );               
                          
    for ( i = 1; i <= nelems; i++ )
    for ( k = 1; k <= nelems; k++ )
    Q[i][k] = double(a[(i-1)*nelems+(k-1)]);
  
 multimat(Q,B1,ht,nelems);
 
 for (i=1,k=1;i<=n;i++)
 {
     if(P[i][i]!=0)
     {
         hst[i]=ht[k];
         k+=1;
     }    
     else 
          hst[i]=0;
 }    
   
  for(i=1;i<=n;i++) 
  {
      if(cells[i][j]!=nodatavalue)
      hstcells[i][j]=hst[i];
	  else
	  hstcells[i][j]=nodatavalue;
  }   
 
 
 for (i=1;i<=nelems;i++)
 free(Q[i]);
 free(Q);
 free(B1);
 free(a);
 free(ipvt);
 free(work);
 
}
 

} //end of J loop 

//Free memory for P, Dx, Dy, hst, B
for(i=1;i<=n;i++)
free(P[i]);
free(P);
free(Dx);
free(Dy);
free(hst);
free(B);










// J direction


//Memory allocation P, Dx, Dy, hst, B with dimension j

Dx = (double *)malloc(((m+1) * sizeof(double)));
Dy =(double *)malloc(((m+1) * sizeof(double)));
hst =(double *)malloc(((m+1) * sizeof(double)));
B = (double *)malloc(((m+1) * sizeof(double )));
P = (double **)malloc(((m+1) * sizeof(double *)));
for(j=1;j<=m;j++)
P[j]=(double *)malloc(((m+1) * sizeof(double)));



//initialize P, Dx, Dy, hst, B with zero







for(i=1;i<=n;i++)
for(j=1;j<=m;j++)
if(cells[i][j]!= nodatavalue)
hwater[i][j]=hstcells[i][j]-cells[i][j];





for(i=1;i<=n;i++)
{

for(j=1;j<=m;j++)
{
   Dx[j]=0; 
   Dy[j]=0; 
   hst[j]=0;
   B[j]=0;
   for(k=1;k<=m;k++)
   P[j][k]=0;
}




for(j=1;j<=m;j++)
{    
   if(i!=1)
   { 
     if((cells[i][j]!= nodatavalue)&&(cells[i-1][j]==nodatavalue))
     K1=0;    
     else if((cells[i][j]!= nodatavalue)&&(cells[i-1][j]!=nodatavalue))
     {
        h1 = 0.5*(hwater[i-1][j]+hwater[i][j]);
		if(j!=m)
			{
			if(cells[i][j+1] == nodatavalue) 
					z1 = hstcells[i][j];
			else
					z1 = hstcells[i][j+1];
        
			if(cells[i-1][j+1] == nodatavalue) 
					z2 = hstcells[i-1][j];
			else
					z2=hstcells[i-1][j+1];
			}
		else
			{
				z1 = hstcells[i][j];
				z2 = hstcells[i-1][j];
			}

		if(j!=1)
			{
			if(cells[i][j-1] == nodatavalue)  
					z3 = hstcells[i][j];
			else
					z3 =hstcells[i][j-1];
			if(cells[i-1][j-1] == nodatavalue) 
					z4 = hstcells[i-1][j];
			else
					z4= hstcells[i-1][j-1];
			}
		else
			{
				z3 = hstcells[i][j];
				z4 = hstcells[i-1][j];	
			}
                
        sqso = (1/(del_x*del_x))*(((hstcells[i][j]-hstcells[i-1][j])*(hstcells[i][j]-hstcells[i-1][j]))+
               (r*r/16)*((z1+ z2- z3- z4)*(z1+ z2- z3- z4)));
     
	  
      if(h1>0)
      K1= (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
      else
      K1=0;
      
	 
	 if(hstcells[i][j]>hstcells[i-1][j])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
	  else if(hstcells[i-1][j]>hstcells[i][j])
		  {
		  if((hwater[i-1][j]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
      else if(hstcells[i-1][j]==hstcells[i][j])
		  K1=0;
	 
	 
	 }
   }
    else
    K1 =0; 
     
     
    if(i!=n)
    { 
     if((cells[i][j]!= nodatavalue)&&(cells[i+1][j]==nodatavalue))
        K2=0;
     else if((cells[i][j]!= nodatavalue)&&(cells[i+1][j]!=nodatavalue))
      {
          h1 = 0.5*(hwater[i][j]+hwater[i+1][j]);
		if(j!=m)
			{
			  if(cells[i][j+1] == nodatavalue)  
						z1= hstcells[i][j];
			  else
						z1=hstcells[i][j+1];
			  if(cells[i+1][j+1] == nodatavalue) 
						z2 = hstcells[i+1][j];
			  else
						z2=hstcells[i+1][j+1];
			}
		else
			{
				z1= hstcells[i][j];
				z2 = hstcells[i+1][j];	
			}

          if(j!=1)
			  {
			  if(cells[i][j-1] == nodatavalue) 
						z3 = hstcells[i][j];
			  else
						z3=hstcells[i][j-1];
			  if(cells[i+1][j-1] == nodatavalue) 
						z4 = hstcells[i+1][j];
			  else
						z4 =hstcells[i+1][j-1];
			  }
          else
			  {
				z3 = hstcells[i][j];
				z4 = hstcells[i+1][j];	
			  }
          
          sqso = (1/(del_x*del_x))*(((hstcells[i+1][j]-hstcells[i][j])*(hstcells[i+1][j]-hstcells[i][j]))+
                 (r*r/16)*((z2+ z1- z4- z3)*(z2+ z1- z4- z3)));
      
      if(h1>0)
      K2= (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
      else
      K2=0;
	  	
 
     if(hstcells[i][j]>hstcells[i+1][j])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
	 else if(hstcells[i+1][j]>hstcells[i][j])
		  {
		  if((hwater[i+1][j]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
    if(hstcells[i+1][j]==hstcells[i][j])
		  K2=0;	
      }                 
    }
    else 
    K2=0;  
       
      
      
      if(cells[i][j] != nodatavalue)
      {
         if(K2==0)
         z1=0;
         else
         z1=K2* (hstcells[i+1][j]-hstcells[i][j]);
         if(K1==0)
         z2=0;
         else
         z2 =K1*(hstcells[i-1][j]-hstcells[i][j]);
         Dx[j] = ((0.5*del_t)/(del_x*del_y))*(z1+z2);
      }
   
  
  
         
               

      
      
      
   if(j!=1)
   { 
      if((cells[i][j]!= nodatavalue)&&(cells[i][j-1]==nodatavalue))
       K1=0;
       
      else if((cells[i][j]!= nodatavalue)&&(cells[i][j-1]!=nodatavalue))
      {
        h1 = 0.5*(hwater[i][j-1]+hwater[i][j]);
		if(i!=n)
			{
			if(cells[i+1][j] == nodatavalue)  
					z1= hstcells[i][j]; 
			else
					z1=hstcells[i+1][j];
			if(cells[i+1][j-1] == nodatavalue) 
					z2 = hstcells[i][j-1];
			else
					z2=hstcells[i+1][j-1];
			}
		else
			{
				z1= hstcells[i][j];
				z2 = hstcells[i][j-1];
			}

		if(i!=1)
			{
			if(cells[i-1][j] == nodatavalue)  
					z3= hstcells[i][j];
			else
					z3=hstcells[i-1][j];
			if(cells[i-1][j-1] == nodatavalue) 
					z4 = hstcells[i][j-1];
			else
					z4=hstcells[i-1][j-1];
			}
		else
			{
				z3= hstcells[i][j];
				z4 = hstcells[i][j-1];
			}
     
        sqso = (1/(del_y*del_y))*(((hstcells[i][j]-hstcells[i][j-1])*(hstcells[i][j]-hstcells[i][j-1]))+
               (r*r/16)*((z1+ z2- z3- z4)*(z1+ z2- z3- z4)));
      
     if(h1>0)
      K1= (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
      else
      K1=0;

	  if(hstcells[i][j]>hstcells[i][j-1])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
	  else if(hstcells[i][j-1]>hstcells[i][j])
		  {
		  if((hwater[i][j-1]<min_h)||(sqso<min_sqso))
				K1=0;
		  else ;
		  }
      else if(hstcells[i][j-1]==hstcells[i][j])
		  K1=0;
	 


      }
   }
   else 
       K1=0;   
      
  
   
  if(j!=m)
  {    
      if((cells[i][j]!= nodatavalue)&&(cells[i][j+1]==nodatavalue))
      K2=0;
      else if((cells[i][j]!= nodatavalue)&&(cells[i][j+1]!=nodatavalue))
      {
        h1 = 0.5*(hwater[i][j+1]+hwater[i][j]);

        if(i!=n)
			{
			if(cells[i+1][j] == nodatavalue)  
				z1= hstcells[i][j];
			else
				z1=hstcells[i+1][j];
			if(cells[i+1][j+1] == nodatavalue) 
				z2= hstcells[i][j+1];
			else
				z2=hstcells[i+1][j+1]; 
			}
		else
			{
				z1= hstcells[i][j];
				z2= hstcells[i][j+1];
			}

		if(i!=1)
			{
			if(cells[i-1][j] == nodatavalue)  
				z3= hstcells[i][j];
			else
				z3=hstcells[i-1][j];
			if(cells[i-1][j+1] == nodatavalue)  
				z4= hstcells[i][j+1];
			else
				z4=hstcells[i-1][j+1];
			}
		else
			{
				z3= hstcells[i][j];
				z4= hstcells[i][j+1];
			}
		
     
        sqso = (1/(del_y*del_y))*(((hstcells[i][j+1]-hstcells[i][j])*(hstcells[i][j+1]-hstcells[i][j]))+
               (r*r/16)*((z2+ z1- z4- z3)*(z2+ z1- z4- z3)));
      
      if(h1>0) 
      K2 = (pow(h1,(5/3.0)))/(nm*pow(sqso,0.25));
	  else
	  K2 =0;
	  
	  
	  if(hstcells[i][j]>hstcells[i][j+1])
		  {
		  if((hwater[i][j]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
	  else if(hstcells[i][j+1]>hstcells[i][j])
		  {
		  if((hwater[i][j+1]<min_h)||(sqso<min_sqso))
				K2=0;
		  else ;
		  }
      else if(hstcells[i][j+1]==hstcells[i][j])
		  K2=0;		


      }
   } 
   else 
   K2=0;     
 
      
        
                     
       
      
      if(cells[i][j] != nodatavalue)
      {
           
         if(K2==0)
         z1=0;
         else
         z1= K2*(hstcells[i][j+1]-hstcells[i][j]);
           
         if(K1==0)
         z2=0;
         else
         z2= K1*(hstcells[i][j-1]-hstcells[i][j]);
         Dy[j] = ((0.5*del_t)/(del_x*del_y))*(z1+z2);
      }

     
     if((j!=1)&&(j!=m)&&(cells[i][j] != nodatavalue))
     {
     P[j][j-1]= -(0.5*del_t*K1)/(del_x*del_y);
     P[j][j] = 1+((0.5*del_t*(K1+K2))/(del_x*del_y));
     P[j][j+1]= -(0.5*del_t*K2)/(del_x*del_y);
     }    
     else if((j==1)&&(j!=m)&&(cells[i][j] != nodatavalue))
     {
     P[j][j] = 1+((0.5*del_t*(K1+K2))/(del_x*del_y)); 
     P[j][j+1]= -(0.5*del_t*K2)/(del_x*del_y);
     }  
     else if((j!=1)&&(j==m)&&(cells[i][j] != nodatavalue)) 
     {
     P[j][j-1]= -(0.5*del_t*K1)/(del_x*del_y);
     P[j][j] = 1+((0.5*del_t*(K1+K2))/(del_x*del_y));
     }
         
  

} //end of J loop






  for(j=1;j<=m;j++) 
  {
      hst[j]=0;
      if(cells[i][j] != nodatavalue)
      {
      if(outlet[i][j]==1)
    
         {
         h1 = hwater[i][j];
         if(h1>min_h)
         K2 = pow(h1,(5/3.0))/nm;
         else K2 = 0;
         sy = (cells[i][j-1]- cells[i][j])/del_y;
         Q2[i][j]= (K2*pow(sy,0.5)*del_x);
         B[j]=Dx[j]+hstcells[i][j]- (K2*pow(sy,0.5)*del_t/(del_y)) +0.5*rain*del_t;
         } 
      else 
        B[j]=Dx[j]+hstcells[i][j]+0.5*rain*del_t; 
      }          
      else
      B[j]= 0;
  }



 

  //Remove zero in diagonals of P Matrix

//count non zero elements in diagonal
nelems = 0;
for (j=1;j<=m;j++)
if(P[j][j]!=0)
nelems = nelems+1;

 if(nelems!=0)

{

//Allocate memory for Q, ht matrix
 Q=(double **)malloc((nelems+1) * sizeof(double *));
 ht = (double *)malloc((nelems+1) * sizeof(double ));
 B1 = (double *)malloc((nelems+1) * sizeof(double ));
 a=(double *)malloc((nelems*nelems) * sizeof(double ));
 ipvt = (int *)malloc((nelems) * sizeof(int));
 work = (double *)malloc((nelems) * sizeof(double));
 for (j=1;j<=nelems;j++)
 Q[j]=(double *)malloc((nelems+1) * sizeof(double ));



 



	  
//initialize  Q matrix
for (j=1;j<=nelems;j++)
{
ht[j]=0;
for (k=1;k<=nelems;k++)
Q[j][k]=0;
}

//Write non zero elements in Q matrix
for (j=1,k=1;j<=m;j++)
{
    if((P[j][j]!=0)&&(j!=1)&&(j!=m))
    {
      Q[k][k]=P[j][j];  
      if(P[j][j-1]!=0) Q[k][k-1]=P[j][j-1];
      if(P[j][j+1]!=0) Q[k][k+1]=P[j][j+1];
      B1[k]=B[j];
      k+=1;
    }
    else if((P[j][j]!=0)&&(j==1)&&(j!=m))
    {
      Q[k][k]=P[j][j];  
      if(P[j][j+1]!=0) Q[k][k+1]=P[j][j+1];
      B1[k]=B[j];
      k+=1;
    }    
    else if((P[j][j]!=0)&&(j!=1)&&(j==m))
    {
      Q[k][k]=P[j][j];  
      if(P[j][j-1]!=0) Q[k][k-1]=P[j][j-1];
      B1[k]=B[j];
      k+=1;
    }    
}    



   
 





 //Inverse Q matrix
 
    for(j=1;j<=nelems;j++)
    for(k=1;k<=nelems;k++)
    a[((j-1))*nelems +(k-1)]= Q[j][k];      
     
    info = dgefa ( a, nelems, nelems, ipvt );

    if ( info != 0 )
    {
       cout << "  Error!  The matrix is nearly singular!\n";
       return 0;
    }        
    job = 11;
    dgedi ( a, nelems, nelems, ipvt, det, work, job );               
                          
    for ( j = 1; j <= nelems; j++ )
    for ( k = 1; k <= nelems; k++ )
    Q[j][k] = double(a[(j-1)*nelems+(k-1)]);
    
    
  
 
 multimat(Q,B1,ht,nelems);
 

 
 for (j=1,k=1;j<=m;j++)
 {
     if(P[j][j]!=0)
     {
         hst[j]=ht[k];
         k+=1;
     }    
     else 
          hst[j]=0;
 }    
   
  for(j=1;j<=m;j++) 
  {
      if(cells[i][j]!=nodatavalue)
      hcells[i][j]=hst[j];
	  else
      hcells[i][j]=nodatavalue;
  }   
  
  

  
  
  //  Free memory for Q B1 ht
  
  
 for (j=1;j<=nelems;j++)
 free(Q[j]);
 free(Q);
 free(B1);
 free(a);
 free(ipvt);
 free(work);  
 }
  


} //end of i loop 



//Free memory for P, Dx, Dy, hst, B
//for(j=1;j<=m;j++)
//free(P[j]);
free(P);
free(Dx);
free(Dy);
free(hst);
free(B);



float k1=0.5;
for(i=1;i<=n;i++)
for(j=1;j<=m;j++)
if(cells[i][j]!= nodatavalue)

    {
        hwater[i][j]=hcells[i][j]-cells[i][j];
  
   
           if(outlet[i][j]!=0)
         
            {
           Qout[i][j] = Q1[i][j]+Q2[i][j];
         
            }

     }
k=0;

for(i=1;i<=n;i++)
for(j=1;j<=m;j++)
{
  
  if(outlet[i][j]!=0)
    {
 // 
    if (i%2!=0)
     {
     
     Qout[i][j]=Qout[i][j]*2;
	 file_open4<<Qout[i][j];
	 
	 //file_open4<<'\t'<<k1;
	 file_open4<<endl;
	 
}
//	 k++;
//	  if(k!=noutlet)
//	  file_open4<<"\t";
//	   else
	 
//	}
     
    
      
    
    
    }       
}
} // End t loop

 

return 0 ;
}// end main







int multimat(double **ptr1,double *ptr2,double *ptr3,int n2)
{
int i,j;
for (i=1;i<=n2;i++)
ptr3[i]=0;
for (i=1;i<=n2;i++)
for (j=1;j<=n2;j++)
ptr3[i]+=ptr1[i][j]*ptr2[j];
return 0;
}     



                      
 /******************************************************
 
 ALL BELOW THIS BELONGS TO INVERSE MATRIX OF LINPACK REFER LINPACK_SUB1.cpp
 
 
 ***********************************************************/
 
 
 int dgefa ( double a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGEFA factors a real general matrix.
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    C++ translation by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].
//    On intput, the matrix to be factored.
//    On output, an upper triangular matrix and the multipliers used to obtain
//    it.  The factorization can be written A=L*U, where L is a product of
//    permutation and unit lower triangular matrices, and U is upper triangular.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int DGEFA, singularity indicator.
//    0, normal value.
//    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
//    but it does indicate that DGESL or DGEDI will divide by zero if called.
//    Use RCOND in DGECO for a reliable indication of singularity.
//
{
  int info;
  int j;
  int k;
  int l;
  double t;
//
//  Gaussian elimination with partial pivoting.
//
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L = pivot index.
//
    l = idamax ( n-k+1, a+(k-1)+(k-1)*lda, 1 ) + k - 1;
    ipvt[k-1] = l;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l-1+(k-1)*lda] == 0.0 )
    {
      info = k;
      continue;
    }
//
//  Interchange if necessary.
//
    if ( l != k )
    {
      t = a[l-1+(k-1)*lda];
      a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda] = t;
    }
//
//  Compute multipliers.
//
    t = -1.0 / a[k-1+(k-1)*lda];

    dscal ( n-k, t, a+k+(k-1)*lda, 1 );
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      t = a[l-1+(j-1)*lda];
      if ( l != k )
      {
        a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = t;
      }
      daxpy ( n-k, t, a+k+(k-1)*lda, 1, a+k+(j-1)*lda, 1 );
    }

  }

  ipvt[n-1] = n;

  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}


void dgedi ( double a[], int lda, int n, int ipvt[], double det[], 
  double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DGEDI computes the determinant and inverse of a matrix factored by DGECO or DGEFA.
//
//  Discussion:
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if DGECO has set 0.0 < RCOND or DGEFA has set INFO == 0.
//
//  Modified:
//
//    17 May 2005
//
//  Author:
//
//    C++ translation by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N], on input, the LU factor information,
//    as output by DGECO or DGEFA.  On output, the inverse
//    matrix if requested.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix A.
//
//    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
//
//    Workspace, double WORK[N].
//
//    Output, double DET[2], the determinant of original matrix if
//    requested.  The determinant = DET[0] * pow ( 10.0, DET[1] )
//    with  1.0 <= abs ( DET[0] ) < 10.0 or DET[0] == 0.0.
//
//    Input, int JOB, specifies what is to be computed.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
{
  int i;
  int j;
  int k;
  int l;
  double t;
//
//  Compute the determinant.
//
  if ( job / 10 != 0 )
  {
    det[0] = 1.0;
    det[1] = 0.0;

    for ( i = 1; i <= n; i++ )
    {
      if ( ipvt[i-1] != i )
      {
        det[0] = -det[0];
      }
      det[0] = det[0] * a[i-1+(i-1)*lda];

      if ( det[0] == 0.0 )
      {
        break;
      }

      while ( r8_abs ( det[0] ) < 1.0 )
      {
        det[0] = det[0] * 10.0;
        det[1] = det[1] - 1.0;
      }
      while ( 10.0 <= r8_abs ( det[0] ) )
      {
        det[0] = det[0] / 10.0;
        det[1] = det[1] + 1.0;
      }
    }
  }
//
//  Compute inverse(U).
//
  if ( ( job % 10 ) != 0 )
  {
    for ( k = 1; k <= n; k++ )
    {
      a[k-1+(k-1)*lda] = 1.0 / a[k-1+(k-1)*lda];
      t = -a[k-1+(k-1)*lda];
      dscal ( k-1, t, a+0+(k-1)*lda, 1 );

      for ( j = k+1; j <= n; j++ )
      {
        t = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = 0.0;
        daxpy ( k, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
      }
    }
//
//  Form inverse(U) * inverse(L).
//
    for ( k = n-1; 1 <= k; k-- )
    {
      for ( i = k+1; i <= n; i++ )
      {
        work[i-1] = a[i-1+(k-1)*lda];
        a[i-1+(k-1)*lda] = 0.0;
      }

      for ( j = k+1; j <= n; j++ )
      {
        t = work[j-1];
        daxpy ( n, t, a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
      }

      l = ipvt[k-1];
      if ( l != k )
      {
        dswap ( n, a+0+(k-1)*lda, 1, a+0+(l-1)*lda, 1 );
      }
    }
  }

  return;
}

//****************************************************************************80

int idamax ( int n, double dx[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    IDAMAX finds the index of the vector element of maximum absolute value.
//
//  Discussion:
//
//    WARNING: This index is a 1-based index, not a 0-based index!
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector to be examined.
//
//    Input, int INCX, the increment between successive entries of SX.
//
//    Output, int IDAMAX, the index of the element of maximum
//    absolute value.
//
{
  double dmax;
  int i;
  int ix;
  int value;

  value = 0;

  if ( n < 1 || incx <= 0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    dmax = r8_abs ( dx[0] );

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < r8_abs ( dx[i] ) )
      {
        value = i + 1;
        dmax = r8_abs ( dx[i] );
      }
    }
  }
  else
  {
    ix = 0;
    dmax = r8_abs ( dx[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < r8_abs ( dx[ix] ) )
      {
        value = i + 1;
        dmax = r8_abs ( dx[ix] );
      }
      ix = ix + incx;
    }
  }

  return value;
}

//****************************************************************************80

void dscal ( int n, double sa, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }

  }

  return;
}

void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539, 
//    ACM Transactions on Mathematical Software, 
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }

  }

  return;
}

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of a R8.
//
//  Modified:
//
//    02 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}

void dswap ( int n, double x[], int incx, double y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DSWAP interchanges two vectors.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539, 
//    ACM Transactions on Mathematical Software, 
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{
  int i;
  int ix;
  int iy;
  int m;
  double temp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    m = n % 3;

    for ( i = 0; i < m; i++ )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }

    for ( i = m; i < n; i = i + 3 )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;

      temp = x[i+1];
      x[i+1] = y[i+1];
      y[i+1] = temp;

      temp = x[i+2];
      x[i+2] = y[i+2];
      y[i+2] = temp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = x[ix];
      x[ix] = y[iy];
      y[iy] = temp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return ;
}

