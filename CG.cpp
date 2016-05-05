//main icindeki onecol ve newcol array leri p-1 den p ye degistirildi.
#include "cg_header.h"


void readInput(char s[]);
void generate_initcols();
void branch(IloNumArrayArray gecmis,IloNumArray dallar);
int heuristic1 (IloNumArray fiyat, IloNumArray ciz);
int heuristic2 (IloNumArray fiyat, IloNumArray ciz);
int node_heuristic1 (IloNumArray fiyat, IloNumArrayArray ek, IloNumArray ciz);
int node_heuristic2 (IloNumArray fiyat, IloNumArrayArray ek, IloNumArray ciz);
int UB(IloNumArrayArray frac_sch, int boyut );
void heuristic1_analysis();
void heuristic2_analysis();

bool integer(float k);


void calculate_ranks();
void calculate_pnumber();
void find_dreg (); // function for assigning regions to districts to be used in heuristic_pp
void region_order();
void find_exactA();

int I,J,T,K;
int d[Imax][Jmax],h[Imax],n[Jmax],r[Jmax],order[Imax][Jmax],rn[Imax],pnumber[Kmax],region_pset[Kmax][Jmax],region_dset[Kmax][Imax],dnumber[Kmax],ordered_regs[Kmax],phar_dset[Jmax][Imax],pd[Jmax];
double avg[Jmax];
char fileName[23],unknown[2];
int numcol,numnodes,numinf,numcgh1,numcgh2,numcg,kol_say,say1,say2,say3,say4,say5,h1imp,h2imp,exactimp,numofcalle,numofcallh1,numofcallh2;
int CoTint,noUsed,UB_imp;
double ToTint;

int init_y[Jmax][Tmax],init_obj[Tmax]; // initial columns
int BestValue,BestSch[Kmax][Tmax],UB_sch[Kmax][Tmax],UBvalue,BestLB;
double together[Jmax][Jmax]; // branching parameters


IloEnv env;
IloNumArrayArray kolon(env,Tmax-1);
IloInt C=Imax;
IloInt P=Jmax;
IloInt G=Tmax;
IloInt B=Kmax;
IloInt Cols=Tmax;
IloNum cputime,pptime,mptime,imp_time,branch_time1,branch_time2,branch_time3,branch_time4,branch_time5,int_time,sch_time;
IloNumArray duty(env, P);
IloNumArray new_pair(env,4);
IloNum cost,heuristic_obj;
IloNumArrayArray sch(env, Kmax);
IloNum sag;
IloTimer timer(env),pptimer(env),mptimer(env),imp_timer(env),branch_timer(env),sch_timer(env),int_timer(env),h1_timer(env),h2_timer(env);
ofstream off;
ofstream outfile,outfile2,summary,results,cpu_mp,col_debug,RMP_sim,UB_sim,single_day;
ofstream pp_analysis;

int i,j,k,t,s,p,bayrak,gear,index,l,flag1,flag2,sayac,ex,rootOpt,Q;
int min_devp,tempi;
int min_dist,min_distj;
double min_dev,tempd;
double min_obj,obj,jchosen;
int flag,ecz1,ecz2,j3,cift;
int size,depth;
double imp_amount_e,imp_amount_h1,imp_amount_h2,av_single_day_e,av_single_day_h1,av_single_day_h2;

FILE *fp;
FILE *fmp_cpu;

int main(int argc, char *argv[])
{
	FILE *finput_multi;
	finput_multi=fopen("instances.txt","r");
	
	//outfile.open("master_report.txt");
	//outfile2.open("sub_report.txt");
	summary.open("summary.txt");
	//results.open("results.txt");
	//off.open("nodes.txt");
	//cpu_mp.open("cpu_analiz_mp.txt");
	//col_debug.open("col_debug.txt");
	//RMP_sim.open("RMP_sim.txt");
	//UB_sim.open("UB_sim.txt");
	//single_day.open("sinle_day.txt");
	pp_analysis.open("pp_exact_results.txt");

	IloNumVarArray2 x(env,C);
	IloNumVarArray y(env,P);
	IloExpr reduced_cost(env,0);
	
	try{

	int ins;
	
	for (ins=1;ins<=InsMax;ins++)
	{
		 // for CPU and RAM analysis
		 timer.restart();
		 cputime=0;
		 pptime=0;
		 mptime=0;
		 imp_time=0;
		 branch_time1=0;
		 branch_time2=0;
		 branch_time3=0;
		 branch_time4=0;
		 branch_time5=0;
		 sch_time=0;
		 int_time=0;
		 numnodes=0;
		 numinf=0;
		 numcg=0;numcgh1=0;numcgh2=0;h1imp=0;h2imp=0;exactimp=0;numofcalle=0;numofcallh1=0;numofcallh2=0;
		 say1=0;
		 say2=0;
		 say3=0;
		 say4=0;
		 say5=0;
		 imp_amount_e=0;imp_amount_h1=0;imp_amount_h2=0;
		 av_single_day_e=0;av_single_day_h1=0;av_single_day_h2=0;
		 CoTint=0;ToTint=0;noUsed=0;UB_imp=0;
		 //---------------------------
		 rootOpt=0;
		 fgets(fileName,24,finput_multi);
		 fgets(unknown,3,finput_multi);
		 readInput(fileName);

		 // auxiliary functions
		 calculate_pnumber();
		 calculate_ranks();
		 readInput(fileName);
		 find_dreg(); // these methods should be opened if heuristic pp is to be used.
		 region_order();
		 find_exactA();
		 //-------------------


		 //off<< "INSTANCE: "<< ins <<endl;
		 //off<<"past kolon size:"<<kolon.getSize()<<endl;
		 // initialize BestValue and BestSch
		 BestValue=M;
		 BestLB=M;
		 for(k=0;k<Kmax;k++)
		 {
			for (t=1;t<Tmax;t++)
				BestSch[k][t]=0;
		 }
		 numnodes=0;//node counter
		 IloNumArray price(env, P); // stores the dual prices of RMP
		 IloNumArray new_col(env, P); // stores the new column found by PP
		 	     		 
		 for(t=0; t < Tmax-1; t++)
		 {
		    kolon[t] = IloNumArray (env, P+1);      
         }

		 for (j=1;j<Jmax;j++) 
		 {
			duty[j]=n[j];
	   	 }

		 // ------------------define the decsision variables-----------------
		 IloNumVarArray lambda(env); 
		 
		 // ------------------build the Master problem model------------------------------
		 IloModel masterPDS(env);
		 IloModel masterIP(env);
		 IloObjective obj_exp = IloAdd(masterPDS, IloMinimize(env));
		 IloRange rng_1 = IloAdd(masterPDS,IloRange(env, Cols-1,IloInfinity ));
		 IloRangeArray rng_2 = IloAdd(masterPDS,IloRangeArray(env, -IloInfinity, duty));
		 
		 //------------------create the initial columns--------------------
		 generate_initcols();
		 for (s = 0; s < Tmax-1; s++) 
		 {
			 kolon[s][0]=init_obj[s+1];
			 kolon[s][1]=1;
			 for (j=2;j<=Jmax;j++) 
				 kolon[s][j]=init_y[j-1][s+1]; 		 
		 }
		 
		 IloNumArray one_col(env,P);
		 for (s = 0; s < Tmax-1; s++) 
		 {
			 for (j=1;j<Jmax;j++) 
				one_col[j]=init_y[j][s+1]; 

			 lambda.add(IloNumVar(obj_exp(init_obj[s+1]) + rng_1(1) + rng_2(one_col) ));
		 }
		 one_col.end();
		 
		 IloCplex MaP(masterPDS);
		 IloCplex MPIP(masterIP);
		 //give the settings
		 MaP.setParam(IloCplex::NumericalEmphasis,1);
		 MaP.setParam(IloCplex::EpAGap,1e-9);
		 MaP.setParam(IloCplex::EpGap,0);
		 MaP.setParam(IloCplex::EpMrk,0.99999);
		 MaP.setParam(IloCplex::EpOpt,1e-9);
		 MaP.setParam(IloCplex::EpPer,1e-8);
		 MaP.setParam(IloCplex::EpRHS,1e-9);
		 MaP.setParam(IloCplex::EpInt,0);

 
		// --------- build the reduced pricing subproblem model------------------------
		
		 for(i=0; i<Imax; i++)
		 {
            x[i] = IloNumVarArray(env, P); 
			for(j=1; j<=rn[i]; j++)
				 x[i][order[i][j]] = IloNumVar(env, 0.0,1.0, ILOFLOAT);
		 }
		 
		 for (j=0;j<Jmax;j++) 
			  y[j]=IloNumVar(env,0,1,ILOBOOL);
		 
		 IloModel subPDS (env);
		 IloObjective ReducedCost = IloAdd(subPDS, IloMinimize(env));
		 
		 for(i=1; i<Imax; i++)
		 {
			for(j=1; j<=rn[i]; j++) 
			{
				IloExpr con1(env);
				con1 = x[i][order[i][j]];
				subPDS.add(con1 <=  y[order[i][j]]);
				con1.end();
			}	
		 }

		 for(i=1; i<Imax; i++)
		 {
			IloExpr con2(env);

			for (j=1;j<=rn[i];j++) 
				con2 += x[i][order[i][j]];
			subPDS.add(con2 >= 1);
			con2.end();				
		 }
	
		for (k=1;k<Kmax;k++)
		{				
			IloExpr con3(env);
			for (j=1;j<Jmax;j++)
			{
				if (r[j]==k)
				   con3 += y[j];
			}			
			subPDS.add(con3 ==  1);
			con3.end();
		}
		
		//-----------------------------------------------------------------------------
		 IloCplex SP(subPDS);
		 SP.setParam(IloCplex::NumericalEmphasis,1);
		 SP.setParam(IloCplex::EpAGap,1e-9);
		 SP.setParam(IloCplex::EpGap,0);
		 SP.setParam(IloCplex::EpMrk,0.99999);
		 SP.setParam(IloCplex::EpOpt,1e-9);
		 SP.setParam(IloCplex::EpPer,1e-8);
		 SP.setParam(IloCplex::EpRHS,1e-9);
		 SP.setParam(IloCplex::EpInt,0);
		 
		//--------------------- SOLUTION PROCEDURE------------------------------
		
	//outfile <<"----------------------------------root node--------------------------------------------"<< endl;
	//outfile2 <<"-------------------------------------root node------------------------------------------- "<< endl;


	//------------------solve the intial model -----------------------
	MaP.exportModel("master.lp");
	mptimer.restart();
	MaP.solve();
	mptime=mptime+mptimer.stop();
	BestLB=MaP.getObjValue();  
	price[0]=MaP.getDual(rng_1);
	for (j=1;j<Jmax;j++)
	{
		 price[j]=MaP.getDual(rng_2[j]);
	}
	//------------------------------------------GEAR-----------------------------------------	
	heuristic1_analysis();
	gear=3;	
	//while(gear<4)
	//{
	  if(gear==1)
	  {
		Q=0;
		for (;;)
	    {

		//-----------------heuristic solution for p-median pricing problem ---------
		numofcallh1++;
		cost=heuristic1(price,new_col);
		heuristic_obj=0;
		heuristic_obj=heuristic_obj+cost-price[0];
		for(j=1;j<Jmax;j++)
			heuristic_obj=heuristic_obj-price[j]*new_col[j];
		 

		//---------------------------report_SP------------------------------
		
		//outfile2 << "Heuristic1 CG" << endl;
	    //outfile2 << "cost: "<<cost<<" pricing obj: "<<heuristic_obj<< endl;

		/*for(j=1;j<Jmax;j++) 
		{
			outfile2 << "y" << j << " = " << new_col[j] << endl;
		}*/
		if (heuristic_obj > -RC_EPS) break;

		numcgh1++;
		av_single_day_h1=av_single_day_h1+cost;
		//single_day<<cost<<" 1"<<endl;

		lambda.add(IloNumVar(obj_exp(cost) + rng_1(1) + rng_2(new_col) ));
		
		index=kolon.getSize();
		kolon.add(IloNumArray (env, P+1));
		kolon[index][0]=cost;
		kolon[index][1]=1;
		for (j=1;j<Jmax;j++)
			 kolon[index][j+1]=new_col[j];	

		MaP.exportModel("master.lp");
		mptimer.restart();
		MaP.solve();
		mptime=mptime+mptimer.stop();

		//------------------get the optimal dual values-----------------------
		 
		price[0]=MaP.getDual(rng_1);
		for (j=1;j<Jmax;j++)
		{
		    price[j]=MaP.getDual(rng_2[j]);
		}

		 //-----------------------------report_MP in heuristic--------------------------------
		 //outfile << "Heuristic1 CG" << endl;
		 //outfile << "Problem status: " << MaP.getStatus() << endl;
		 //outfile << "Objective of the RMP: " << MaP.getObjValue() << endl;
		 //outfile << "Values of decsision variables:" << endl;

		 /*for ( s = 0; s < lambda.getSize(); s++) {
			outfile << "lambda" << s << " = " << MaP.getValue(lambda[s]) << endl;
		 }
		 outfile << "Prices:" <<endl;

		 for (j = 0; j < Jmax; j++) {
			  outfile << "con_" << j << " = " << price[j] << endl;
		 }*/

		if(MaP.getObjValue()<BestLB)
		{
			 imp_amount_h1=imp_amount_h1+(BestLB-MaP.getObjValue());
			 h1imp++;
			 BestLB=MaP.getObjValue(); 
			 //RMP_sim<<BestLB<<" 1"<<endl;
			 Q=0; 
		}
		else //break;
		  Q++;

		if(Q>25) break;
		
	    }// heuristic1 end
		gear++;
	  }	//gear==1 end
		//----------------------- heuristic solution 2 to pricing problem---------
	 
	 if(gear==2)
	 {
	    bayrak=0;
		Q=0;
		for (;;)
	    {
		
		//-----------------heuristic solution for p-median pricing problem ---------
		numofcallh2++;
		cost=heuristic2(price,new_col);
		heuristic_obj=0;
		heuristic_obj=heuristic_obj+cost-price[0];
		for(j=1;j<Jmax;j++)
			heuristic_obj=heuristic_obj-price[j]*new_col[j];
		 
		 

		//---------------------------report_SP------------------------------
		
		//outfile2 << "Heuristic2 CG" << endl;
	    //outfile2 << "cost: "<<cost<<" pricing obj: "<<heuristic_obj<< endl;

		/*for(j=1;j<Jmax;j++) 
		{
			outfile2 << "y" << j << " = " << new_col[j] << endl;
		}*/
		if (heuristic_obj > -RC_EPS) break;

		numcgh2++;
		av_single_day_h2=av_single_day_h2+cost;
		//single_day<<cost<<" 2"<<endl;

		lambda.add(IloNumVar(obj_exp(cost) + rng_1(1) + rng_2(new_col) ));
		
		index=kolon.getSize();
		kolon.add(IloNumArray (env, P+1));
		kolon[index][0]=cost;
		kolon[index][1]=1;
		for (j=1;j<Jmax;j++)
			 kolon[index][j+1]=new_col[j];		 
		
		MaP.exportModel("master.lp");
		mptimer.restart();
		MaP.solve();
		mptime=mptime+mptimer.stop();

		//------------------get the optimal dual values-----------------------
		 
		 price[0]=MaP.getDual(rng_1);
		 for (j=1;j<Jmax;j++)
		 {
			 price[j]=MaP.getDual(rng_2[j]);
		 }

		 //-----------------------------report_MP in heuristic--------------------------------
		 //outfile << "Heuristic2 CG" << endl;
		 //outfile << "Problem status: " << MaP.getStatus() << endl;
		 //outfile << "Objective of the RMP: " << MaP.getObjValue() << endl;
		 //outfile << "Values of decsision variables:" << endl;

		 /*for ( s = 0; s < lambda.getSize(); s++) {
			outfile << "lambda" << s << " = " << MaP.getValue(lambda[s]) << endl;
		 }
		 outfile << "Prices:" <<endl;

		 for (j = 0; j < Jmax; j++) {
			  outfile << "con_" << j << " = " << price[j] << endl;
		 }*/

		 if(MaP.getObjValue()<BestLB)
		 {
		  
			 imp_amount_h2=imp_amount_h2+(BestLB-MaP.getObjValue());
			 h2imp++;  
			 BestLB=MaP.getObjValue();
			 gear=1;
			 bayrak=1;
			 //RMP_sim<<BestLB<<" 2"<<endl;
			 break;
		 }
		 
		 else //break;
			Q++;

		 if (Q>25) break;
		 

		
	    }// heuristic2 end 
		if(bayrak==1)
			gear=1;
		else
			gear++;
	 }	//gear2 end
		//-----------------exact CG iteration-----------------------------------------------------
	 kol_say=0;
     if(gear==3)
	 {
		bayrak=0;
		for (;;)
	    { 
		 
		 //----------------------------------------------------------------------
		 reduced_cost.clear();
		 for( j=1; j< Jmax; j++)
		 { 		         				
			reduced_cost-=price[j]*y[j];
		 }
		 // for the reduced problem------------------
		 for( i=1; i< Imax; i++)
		 { 
			 for (j=1;j<=rn[i];j++)
		 		reduced_cost += x[i][order[i][j]]*d[i][order[i][j]] * h[i];		         				
		 }
		 //------------------------------------------
		 reduced_cost-=price[0];
		 ReducedCost.setExpr(reduced_cost);  
		 
		 numofcalle++;
		 //SP.exportModel("pricing.lp");
		 pptimer.restart();
		 SP.solve();
		 pptime=pptime+pptimer.stop();

		 // PP_analysis output
		 pp_analysis << numofcalle << " " << price[0] << " ";
		 for(j=1;j<Jmax;j++) pp_analysis << price[j] << " ";
		 pp_analysis << pptimer.stop() << " " << SP.getObjValue() << endl;

		// -------------get the new column------------------------
		for(j=1;j<Jmax;j++) 
		{
			if (SP.getValue(y[j]) < RC_EPS)
				new_col[j] = 0;
			else
				new_col[j] = 1;
		}
		
		cost=0;

		// for reduced pricing-----------------------
		for (i=1;i<Imax;i++)
		{
			min_dist=M;
			for (j=1;j<=rn[i];j++)
			{
				if (new_col[order[i][j]]==1 && d[i][order[i][j]]< min_dist)
				{
					min_distj=order[i][j];
					min_dist=d[i][order[i][j]];
				}
			}
			cost=cost+d[i][min_distj]*h[i];
		}
		//--------------------
		
		//---------------------------report_SP------------------------------
		
		/*outfile2 << "Problem status: " << SP.getStatus() << endl;
		outfile2 << "Objective of the PP: " << SP.getObjValue() << endl;
	    outfile2 << "New column:" << endl;
		for(j=1;j<Jmax;j++) 
		{
			outfile2 << "y" << j << " = " << new_col[j] << endl;
		}*/

		if (SP.getValue(ReducedCost) > -RC_EPS) break;

		numcg++;
		av_single_day_e=av_single_day_e+cost;
		//single_day<<cost<<" 3"<<endl;

		//outfile2 << cost << endl;
		lambda.add(IloNumVar(obj_exp(cost) + rng_1(1) + rng_2(new_col) ));
		kol_say++;

		index=kolon.getSize();
		kolon.add(IloNumArray (env, P+1));
		kolon[index][0]=cost;
		kolon[index][1]=1;
		for (j=1;j<Jmax;j++)
			 kolon[index][j+1]=new_col[j];

		MaP.exportModel("master.lp");
		mptimer.restart();
		MaP.solve();
		mptime=mptime+mptimer.stop(); 
		//cpu_mp<<mptimer.stop()<<" "<<mptime<<endl;
		 
		 //------------------get the optimal dual values-----------------------
		 
		 price[0]=MaP.getDual(rng_1);
		 for (j=1;j<Jmax;j++)
		 {
			 price[j]=MaP.getDual(rng_2[j]);
		 }
		 //-----------------------------report_MP--------------------------------
	
		 /*outfile << "Problem status: " << MaP.getStatus() << endl;
		 outfile << "Objective of the RMP: " << MaP.getObjValue() << endl;
		 outfile << "Values of decsision variables:" << endl;

		 for ( s = 0; s < lambda.getSize(); s++) {
			outfile << "lambda" << s << " = " << MaP.getValue(lambda[s]) << endl;
		 }
		 
		 outfile << "Prices:" <<endl;

		 for (j = 0; j < Jmax; j++) {
			  outfile << "con_" << j << " = " << price[j] << endl;
		 }*/

		 if(MaP.getObjValue()<BestLB)
		 {
			  imp_amount_e=imp_amount_e+(BestLB-MaP.getObjValue());
			  exactimp++;
			  BestLB=MaP.getObjValue(); 
			  //gear=1;
			  //bayrak=1;
			  //RMP_sim<<BestLB<<" 3"<<endl;
			  //break;	 
		 }

		

	  }	// CG iteration
	  if(bayrak==1 ) 
	  		gear=1;
	  else
		  gear=4;
	
	}//gear3 end
//} //while loop end

	col_debug<<kol_say<<endl;
		
	// --------------write the optimal schedule in the root node----------------------------
	/*off<<"root kolons"<<endl;
	for (s=0;s<kolon.getSize();s++)
	{
		for (j=0;j<=Jmax;j++)
			off<<kolon[s][j]<<" ";
		 off<<endl;
	}*/
	IloNum count; // keeps number of basic variables

	count=0;
	for (s=0;s<lambda.getSize();s++)
	{
		if (MaP.getValue(lambda[s]) > 0)
		{
			count=count+1;
		}		
	}
	
	//IloNumArrayArray sch(env, B); 		    		 
	for(k=0; k < Kmax; k++)
    {
		sch[k] = IloNumArray (env, 5000);      
    }
	
	count=0;
	for (s=0;s<lambda.getSize();s++)
	{
	  if (MaP.getValue(lambda[s])>0)	
	  {
		sch[0][count]=MaP.getValue(lambda[s]);

		for (j=1;j<Jmax;j++)
		{	
			if (kolon[s][j+1]==1) 
			   sch[r[j]][count]= j;
		}
		count=count+1;
	  }
	}	

	cputime=timer.stop();
	timer.start();
	/*results <<"Root node optimal RMP for ins:" << ins<<"  Optimal solution value: "<<MaP.getObjValue()<<endl;
	for (s=0;s<count;s++)
	{
			results << sch[0][s]<<"	";
	}
	results << endl;
	
	for (k=1;k<Kmax;k++)
	{
		for (s=0;s<count;s++)	
			 results << sch[k][s] <<"	";
		results << endl;
	}*/
	summary<<"ROOT NODE STATISTICS"<<endl;
	summary<<"ins"<<ins<<" :"<<endl;
	summary<<"Root node cost: "<<MaP.getObjValue()<<endl;
	summary<<"Root node cpu: "<<cputime<<endl;
	summary<<"PP time in the root node:"<<pptime<<endl;
	summary <<"Number of calls to exact pp: "<<numofcalle<<endl;
	summary <<"Number of calls to h1 pp: "<<numofcallh1<<endl;
	summary <<"Number of calls to h2 pp: "<<numofcallh2<<endl;
	summary<<"Number of columns in the root node: "<<lambda.getSize()<<endl;
	summary <<"Number of columns by exact pp: "<<numcg<<endl;
	summary <<"Number of columns by h1 pp: "<<numcgh1<<endl;
	summary <<"Number of columns by h2 pp: "<<numcgh2<<endl;
	summary <<"# of exact imp: "<<exactimp<<endl;
	summary <<"# of h1 imp: "<<h1imp<<endl;
	summary <<"# of h2 imp: "<<h2imp<<endl;
	summary <<"average RMP improvement by exact pp: "<<(double)imp_amount_e/exactimp<<endl;
	summary <<"average sinle day cost found by exact pp: "<<(double)av_single_day_e/numcg<<endl;
	summary <<"average RMP improvement by h1 pp: "<<(double)imp_amount_h1/h1imp<<endl;
	summary <<"average sinle day cost found by h1 pp: "<<(double)av_single_day_h1/numcgh1<<endl;
	summary <<"average RMP improvement by h2 pp: "<<(double)imp_amount_h2/h2imp<<endl;
	summary <<"average sinle day cost found by h2 pp: "<<(double)av_single_day_h2/numcgh2<<endl;


	flag=0;	
	for (s=0;s<count;s++)	
	{
		if (integer(sch[0][s]) == 0)
			flag=1;
	}
	if (flag==0)
		summary<<"Optimal integer at root node due LP:) "<<endl;
	/*else
	{
	fp = fopen("master.lp","r+");
	fseek(fp,-5, SEEK_END); 
	fprintf(fp,"General\n");
	for (s=0;s<=lambda.getSize();s++)
		fprintf(fp,"x%d ",s+1);
	fprintf(fp,"\n");
	fprintf(fp,"END\n");
	fclose(fp); 

	MPIP.importModel(masterIP, "master.lp");
    MPIP.extract(masterIP);
	MPIP.solve();
	cputime=timer.stop();
	timer.start();
	summary<<"Root node incumbent: "<<MPIP.getObjValue()<<" found at time "<<cputime<<endl;
	BestValue=MPIP.getObjValue();
	if (MPIP.getObjValue()== MaP.getObjValue())
	{
		summary<<"Optimal integer at root node due IP:) "<<endl;
		rootOpt=1;
	}
	//UB_sim<<BestValue<<endl;
	}*/
	//-------------------------------------------------------branching time------------------------------------------------------------
	//single_day<<"branching part"<<endl;


	IloNumArrayArray used(env, 1); // keeps the history of branched pairs
    for(i=0; i < 1; i++)
	{
		used[i] = IloNumArray (env, 4);      
	}	

	/*flag=0;	
	for (s=0;s<count;s++)	
	{
		if (integer(sch[0][s]) == 0)
			flag=1;
	}*/	
	//results <<"flag is: "<<flag<<endl;
	
	for(j=1;j<Jmax;j++)
	{
		for(p=1;p<Jmax;p++)
			together[j][p]=0;
	}
	for (s=0;s<count;s++)
	{
		for(p=1;p<Kmax-1;p++)
		{
			ecz1=sch[p][s];
			for (j=p+1;j<Kmax;j++)
			{
				ecz2=sch[j][s];
				if (ecz2<ecz1)
					together[ecz2][ecz1]=together[ecz2][ecz1]+sch[0][s];
				else
					together[ecz1][ecz2]=together[ecz1][ecz2]+sch[0][s];
			}
		}		
	}
	
	/*results<<"together values"<<endl;
	for(j=1;j<Jmax;j++)
	{
		for(p=1;p<Jmax;p++)
		{
			results<< together[j][p] <<"	";
		}
		results<<endl;
	}*/

	// make branching only if flag=1
	if (flag==1 && rootOpt==0)
	{
	 			
	 // lets make the branching
	 IloNumArray brvar(env,3);
	 for(j=1;j<Jmax;j++)
	 {
		for(p=j+1;p<Jmax;p++)
		{
			if(integer(together[j][p]) == 0)
			{
				brvar[0]=j;
				brvar[1]=p;
				brvar[2]=together[j][p];
			}
		}
	 }
	 //branch(used,brvar);
	 brvar.end();

	 /*results<<"Optimal sol: ";
	 results<<BestValue<<endl;
	 for (s=0;s<numcol;s++)	
			 results<<BestSch[0][s]<<"	 ";
		results<<endl;
	 for (k=1;k<Kmax;k++)
	 {
	 	for (s=0;s<numcol;s++)	
			 results<<BestSch[k][s]<<"	 ";
		results<<endl;
	 }*/

	cputime=timer.stop();
	summary<<"OVERALL STATISTICS"<<endl;
	summary <<"Best Obj Value: "<< BestValue <<endl ;
	summary <<"Total Cpu Time: "<< cputime <<endl;
	summary <<"Total PP Time: "<<pptime <<endl;
	summary <<"Total MP Time: "<<mptime <<endl;
	summary <<"Total imp Time: "<<imp_time <<endl;
	summary <<"sch Time: "<<sch_time <<endl;
	summary <<"int Time: "<<int_time <<endl;
	summary <<"branch_time1: "<<branch_time1 <<endl;
	summary <<"branch_time2: "<<branch_time2 <<endl;
	summary <<"branch_time3: "<<branch_time3 <<endl;
	summary <<"branch_time4: "<<branch_time4 <<endl;
	summary <<"branch_time5: "<<branch_time5 <<endl;
	summary <<"Number of nodes: "<<numnodes << endl;
	summary <<"Number of calls to exact pp: "<<numofcalle<<endl;
	summary <<"Number of calls to h1 pp: "<<numofcallh1<<endl;
	summary <<"Number of calls to h2 pp: "<<numofcallh2<<endl;
	summary <<"Total Number of columns: "<<kolon.getSize()<<endl;
	summary <<"Number of columns by exact pp: "<<numcg<<endl;
	summary <<"Number of columns by h1 pp: "<<numcgh1<<endl;
	summary <<"Number of columns by h2 pp: "<<numcgh2<<endl;
	summary <<"Number of infeasible nodes: "<<numinf<<endl; 
	summary <<"Number of calls to UB: "<<noUsed<<endl;
	summary <<"Number of improvements by UB: "<<UB_imp<<endl;
	summary <<"Average # of integer kolons if not all columns are integer: "<<(double)(ToTint/CoTint)<<endl;
	summary <<"say1: "<<say1<<endl;
	summary <<"say2: "<<say2<<endl;
	summary <<"say3: "<<say3<<endl;
	summary <<"say4: "<<say4<<endl;
	summary <<"say5: "<<say5<<endl;
	summary <<"# of exact imp: "<<exactimp<<endl;
	summary <<"# of h1 imp: "<<h1imp<<endl;
	summary <<"# of h2 imp: "<<h2imp<<endl;
	summary <<"average RMP improvement by exact pp: "<<(double)imp_amount_e/exactimp<<endl;
	summary <<"average sinle day cost found by exact pp: "<<(double)av_single_day_e/numcg<<endl;
	summary <<"average RMP improvement by h1 pp: "<<(double)imp_amount_h1/h1imp<<endl;
	summary <<"average sinle day cost found by h1 pp: "<<(double)av_single_day_h1/numcgh1<<endl;
	summary <<"average RMP improvement by h2 pp: "<<(double)imp_amount_h2/h2imp<<endl;
	summary <<"average sinle day cost found by h2 pp: "<<(double)av_single_day_h2/numcgh2<<endl;

	}//end of if statement
	else
	{
		BestValue=MaP.getObjValue();
		for (s=0;s<count;s++)	
				 BestSch[0][s]=sch[0][s];
		for (k=1;k<Kmax;k++)
		{
			for (s=0;s<count;s++)	
				 BestSch[k][s]=sch[k][s];
		}
	}

	//---------------------------------------------------------------------------------------------------------------------------------
	
	} // instance iteration
	env.end();
	
 } catch (IloException e){cout <<e;}
 
 fclose (finput_multi);
 //outfile.close();
 //outfile2.close();
 summary.close();
 //results.close();
 //off.close();
 //cpu_mp.close();
 //RMP_sim.close();
 //UB_sim.close();
 //single_day.close();
 pp_analysis.close();

 return 0;
}

void calculate_ranks()
{
    //printf("calculate_ranks\n");
    int j,mindist_j;
    int p[Kmax];
    //int temp_d[Imax][Jmax];
	int mindist;


    /*for(i=1;i<Imax;i++)
	{
        for (j=1;j<Jmax;j++)
             temp_d[i][j]=d[i][j];
	}*/

    for(k=0;k<Kmax;k++)
    {
         p[k]=0;
    }
     //fprintf(fout,"Ranks:\n");
	rn[0]=0;
     for (i=1;i<Imax;i++)
     {
         rn[i]=1;
		 
         while(1)
         {
                 mindist=M;
				 
                 for (j=1;j<Jmax;j++)
                 {
                     if(d[i][j] < mindist)//---------------------------------------
                     {
						mindist=d[i][j];//----------------------
						mindist_j=j;
                     }
                 }

                 order[i][rn[i]]=mindist_j;
                 rn[i]=rn[i]+1;
                 d[i][mindist_j]=M;//------------------------
                 k=r[mindist_j];
                 p[k]=p[k]+1;
                 if(p[k]== pnumber[k] || rn[i]== Jmax )
					break;
         }

         rn[i]=rn[i]-1;
         for(k=0;k<Kmax;k++)
         {
         p[k]=0;
         }
		 
     }
	 /*off<<"orders"<<endl;
	 for (i=1;i<Imax;i++)
	 {
		for(j=1;j<=rn[i];j++)
			off<<order[i][j]<<" ";
		off<<endl;
	 }*/
     
}
void calculate_pnumber()
{
    //printf("calculate_pnumber\n");
    
	for (k=1; k<Kmax ; k++)
    {
        pnumber[k]=0;
        for (j=1; j<Jmax; j++)
        {
            if (r[j]==k)
            {
                pnumber[k]=pnumber[k]+1;
                region_pset[k][pnumber[k]]=j;
            }
        }
    }

}

void find_dreg ()
{
	
	int dreg[Imax];// regions for districts

	// find the regions of each district
	for(i=1;i<Imax;i++)
	{
		dreg[i]=r[order[i][1]];
	}
	for (k=1; k<Kmax ; k++)
    {
        dnumber[k]=0;
        for (i=1; i<Imax; i++)
        {
            if (dreg[i]==k)
            {
                dnumber[k]=dnumber[k]+1;
                region_dset[k][dnumber[k]]=i;
            }
        }
    }

}

void branch(IloNumArrayArray gecmis,IloNumArray dallar)
{
	int cg;
	int isBetter,isFractional;
	

for (cg=1;cg<=2;cg++)
{
	
	IloEnv env2;
	IloNumVarArray2 x(env2,C);
	IloNumVarArray y(env2,P);
	
	for(i=0; i<Imax; i++)
		 {
            x[i] = IloNumVarArray(env2, P); 
			for(j=1; j<=rn[i]; j++)
				 x[i][order[i][j]] = IloNumVar(env2, 0.0,1.0, ILOFLOAT);
		 }
		 
		 for (j=0;j<Jmax;j++) 
			  y[j]=IloNumVar(env2,0,1,ILOBOOL);



	IloExpr reduced_cost(env2,0);
	
	
	if (cputime > 9000)
		break;
	numnodes++;
	//RMP_sim<<"node "<<numnodes<<"-----------------------------"<<endl;
	BestLB=M;
	
	size=kolon.getSize();
	depth=Jmax+gecmis.getSize();
	

	IloNumArray dal(env2,3);
	for(i=0;i<=2;i++)
		dal[i]=dallar[i];

	IloNumArrayArray hist(env,gecmis.getSize()); // copy the used pairs until now
	for(i=0; i < gecmis.getSize(); i++)
	{
		hist[i] = IloNumArray (env, 4);      
	}

	for(s=0;s<gecmis.getSize();s++)
	{
		for (j=0;j<4;j++)
			hist[s][j]=gecmis[s][j];
	}

	IloModel master1(env2); 
	IloObjective obj_exp = IloAdd(master1, IloMinimize(env2));
	IloNumVarArray decs(env2); 
	IloModel sub1(env2); 
	IloObjective obj1 = IloAdd(sub1, IloMinimize(env2));
	
	IloRangeArray rng_3(env2,depth); 
	IloNumArray one_col(env2,depth);
	IloNumArray price(env2, depth);
	IloNumArray new_col(env2,depth);
	
	
	//-----------------------------------potential improvement area----------------------------------------
	imp_timer.restart();

	rng_3[0]=IloAdd(master1,IloRange(env2, Tmax-1, IloInfinity));
	for (j=1;j<Jmax;j++)
	    rng_3[j] = IloAdd(master1,IloRange(env2, -IloInfinity, duty[j]));
	
	if((depth-Jmax)>1)
	{
	 for (j=Jmax;j<depth-1;j++)
	 {
		if(hist[j-Jmax+1][0]==1)
		 rng_3[j]=IloAdd(master1,IloRange(env2, -IloInfinity, hist[j-Jmax+1][3] ));

		if(hist[j-Jmax+1][0]==(-1))
		 rng_3[j]=IloAdd(master1,IloRange(env2, hist[j-Jmax+1][3], IloInfinity )); 
	 }
	}

	if (cg==1)
	{
	sag=floor(dal[2]);
	rng_3[depth-1]=IloAdd(master1,IloRange(env2, -IloInfinity, sag ));
	}
	else
	{
	sag=ceil(dal[2]);
	rng_3[depth-1]=IloAdd(master1,IloRange (env2, sag, IloInfinity));
	}


	for (s = 0; s < size; s++) 
	 {
		one_col[0]=1;	
		for (j=1;j<J+1;j++) 
		{
		    one_col[j]=kolon[s][j+1]; 
		}
		if((depth-Jmax)>1)
		{
			for (j=1;j<hist.getSize();j++)
			{
				if(kolon[s][hist[j][1]+1]==1 && kolon[s][hist[j][2]+1]==1)
				{
					one_col[J+j]=1;	
				}
				else
					one_col[J+j]=0;
			}
		}

		flag1=0;
		flag2=0;
		
		if (kolon[s][dal[0]+1]==1)
			flag1=1;
		if (kolon[s][dal[1]+1]==1)
			flag2=1;
			
		if (flag1==1 && flag2==1)
		{
			one_col[depth-1]=1;
		}
		else
		{
			one_col[depth-1]=0;
		}
		
	    decs.add(IloNumVar(obj_exp(kolon[s][0]) + rng_3(one_col) ));
		
	 }
	imp_time=imp_time+imp_timer.stop();
	//-----------------------end of improvement area----------------------------------
	//---------------reduced pricing model-------------------------
	for(i=1; i<Imax; i++)
		 {
			for(j=1; j<=rn[i]; j++) 
			{
				IloExpr con1(env2);
				con1 = x[i][order[i][j]];
				sub1.add(con1 <=  y[order[i][j]]);
				con1.end();
			}	
		 }

		 for(i=1; i<Imax; i++)
		 {
			IloExpr con2(env2);

			for (j=1;j<=rn[i];j++) 
				con2 += x[i][order[i][j]];
			sub1.add(con2 >= 1);
			con2.end();				
		 }
	
		for (k=1;k<Kmax;k++)
		{				
			IloExpr con3(env2);
			for (j=1;j<Jmax;j++)
			{
				if (r[j]==k)
				   con3 += y[j];
			}			
			sub1.add(con3 ==  1);
			con3.end();
		}
		//---------------------------------------------------
	if(cg==1)
	new_pair[0]=1;// 1 means less than, -1 means greater than
	else
	new_pair[0]=-1;// 1 means less than, -1 means greater than

	new_pair[1]=dal[0];
	new_pair[2]=dal[1];
	if(cg==1)
		new_pair[3]=floor(dal[2]);
	else
		new_pair[3]=ceil(dal[2]);
	hist.add(new_pair);
	
	
	IloCplex MP(master1);
	MP.setParam(IloCplex::NumericalEmphasis,1);
	MP.setParam(IloCplex::EpAGap,1e-9);
	MP.setParam(IloCplex::EpGap,0);
	MP.setParam(IloCplex::EpMrk,0.99999);
	MP.setParam(IloCplex::EpOpt,1e-9);
	MP.setParam(IloCplex::EpPer,1e-8);
	MP.setParam(IloCplex::EpRHS,1e-9);
	MP.setParam(IloCplex::EpInt,0);

	IloCplex SP(sub1);
	SP.setParam(IloCplex::NumericalEmphasis,1);
	SP.setParam(IloCplex::EpAGap,1e-9);
	SP.setParam(IloCplex::EpGap,0);
	SP.setParam(IloCplex::EpMrk,0.99999);
	SP.setParam(IloCplex::EpOpt,1e-9);
	SP.setParam(IloCplex::EpPer,1e-8);
	SP.setParam(IloCplex::EpRHS,1e-9);
	SP.setParam(IloCplex::EpInt,0);

	//----------------------solve the initial model------------------------------------
	mptimer.restart();
	MP.solve();
	mptime=mptime+mptimer.stop();
	 
	if (MP.getStatus()== MP.Infeasible)
	{
		numinf++;
	}
		//------------------get the optimal dual values-------
	if(MP.getStatus()!= MP.Infeasible)
	{
		price[0]=MP.getDual(rng_3[0]);
		
		for (j=1;j<Jmax;j++)
		{ 
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}

		for (j=Jmax;j<depth;j++)
		{
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}
	}
	

	//---------------------------------------------------------------------------------------------GEAR-----------------------------------
	if(MP.getStatus()!= MP.Infeasible)
	{
	BestLB=MP.getObjValue();
	gear=3;	
	//while(gear<4)
	//{
	  if(gear==1)
	  {
		Q=0;	 
		for (;;)
		{
		
	//-----------------heuristic solution for p-median pricing problem ---------
	numofcallh1++;
	cost=node_heuristic1(price,hist,new_col);
	heuristic_obj=0;
	heuristic_obj=heuristic_obj+cost-price[0];

	for(j=1;j<Jmax;j++)
	{
		heuristic_obj=heuristic_obj-price[j]*new_col[j];
		
		for (l=1;l<hist.getSize();l++)
			{
				
				if(hist[l][1]==j || hist[l][2]==j)
				{
					if(hist[l][0]==1)
						heuristic_obj=heuristic_obj-price[Jmax+l-1];
					else
						heuristic_obj=heuristic_obj+price[Jmax+l-1];
				}  
			}
		
	 }

	 for (j=Jmax;j<depth;j++) 
		{
			flag1=0;flag2=0;
			if (new_col[hist[j-Jmax+1][1]] ==1 )
				flag1=1;
			if (new_col[hist[j-Jmax+1][2]] ==1 )
				flag2=1;
		   
		    if(flag1==1 && flag2==1)
			   new_col[j]=1;
			else
			   new_col[j]=0;
		
		}
	 //---------------------------report_SP------------------------------
		/*outfile2 << "Heuristic1 CG" << endl;
		outfile2 << "cost: "<<cost<<" Objective of the PP: " << heuristic_obj << endl;
	    outfile2 << "New column:" << endl;
		for(j=1;j<Jmax;j++) 
		{
			outfile2 << "y" << j << " = " << new_col[j] << endl;
		}*/

		
		if ( heuristic_obj > -RC_EPS) break;

		numcgh1++;
		av_single_day_h1=av_single_day_h1+cost;
		//single_day<<cost<<" 1"<<endl;
		decs.add(IloNumVar(obj_exp(cost) + rng_3(new_col) ));

	    index=kolon.getSize();
		kolon.add(IloNumArray (env, P+1)); 
		kolon[index][0]=cost;
		kolon[index][1]=1;
		for (j=1;j<Jmax;j++)
		    kolon[index][j+1]=new_col[j];

		mptimer.restart();
		MP.solve();
		mptime=mptime+mptimer.stop();
		 
		if (MP.getStatus()== MP.Infeasible)
		{
			numinf++;
			break;
		}

		//------------------get the optimal dual values-------
		 
		price[0]=MP.getDual(rng_3[0]);
		
		for (j=1;j<Jmax;j++)
		{ 
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}

		for (j=Jmax;j<depth;j++)
		{
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}
		 /*outfile << "Heuristic1 CG" << endl;
		 outfile << "Problem status: " << MP.getStatus() << endl;
		 outfile << "Objective of the RMP: " << MP.getObjValue() << endl;
		 outfile << "Values of decsision variables:" << endl;

		 for ( s = 0; s < decs.getSize(); s++) {
			 outfile << "lambda" << s << " = " << MP.getValue(decs[s]) <<" reduced cost = "<<MP.getReducedCost(decs[s])<< endl;
		 }
		 
		 outfile << "Prices:" <<endl;

		 for (j = 0; j < depth; j++) {
			  outfile << "con_" << j << " = " << price[j]<< endl; 
		 }*/

		 if(MP.getObjValue()<BestLB)
		 {
			imp_amount_h1=imp_amount_h1+(BestLB-MP.getObjValue());
			h1imp++;
			BestLB=MP.getObjValue(); 
			//RMP_sim<<BestLB<<" 1"<<endl;
			Q=0;
		 }
		else //break;
			Q++;

		if(Q>25) break;
		 

	}// end of heuristic pp1
	    
		gear++;
	}	//gear==1 end
	 
	 if(gear==2)
	 {
	    bayrak=0;
		Q=0;
		for (;;)
		 {
 
	//-----------------heuristic solution for p-median pricing problem ---------
	numofcallh2++;
	cost=node_heuristic2(price,hist,new_col);
	heuristic_obj=0;
	heuristic_obj=heuristic_obj+cost-price[0];

	for(j=1;j<Jmax;j++)
	{
		heuristic_obj=heuristic_obj-price[j]*new_col[j];
		
		for (l=1;l<hist.getSize();l++)
			{
				
				if(hist[l][1]==j || hist[l][2]==j)
				{
					if(hist[l][0]==1)
						heuristic_obj=heuristic_obj-price[Jmax+l-1];
					else
						heuristic_obj=heuristic_obj+price[Jmax+l-1];
				}  
			}
		
	 }

	 for (j=Jmax;j<depth;j++) 
		{
			flag1=0;flag2=0;
			if (new_col[hist[j-Jmax+1][1]] ==1 )
				flag1=1;
			if (new_col[hist[j-Jmax+1][2]] ==1 )
				flag2=1;
		   
		    if(flag1==1 && flag2==1)
			   new_col[j]=1;
			else
			   new_col[j]=0;
		
		}
	 //---------------------------report_SP------------------------------
		/*outfile2 << "Heuristic2 CG" << endl;
		outfile2 << "cost: "<<cost<<" Objective of the PP: " << heuristic_obj << endl;
	    outfile2 << "New column:" << endl;
		for(j=1;j<Jmax;j++) 
		{
			outfile2 << "y" << j << " = " << new_col[j] << endl;
		}*/

		
		if ( heuristic_obj > -RC_EPS) break;
		
		numcgh2++;
		av_single_day_h2=av_single_day_h2+cost;
		//single_day<<cost<<" 2"<<endl;
		decs.add(IloNumVar(obj_exp(cost) + rng_3(new_col) ));

	    index=kolon.getSize();
		kolon.add(IloNumArray (env, P+1)); 
		kolon[index][0]=cost;
		kolon[index][1]=1;
		for (j=1;j<Jmax;j++)
		    kolon[index][j+1]=new_col[j];

		mptimer.restart();
		MP.solve();
		mptime=mptime+mptimer.stop();
		 
		if (MP.getStatus()== MP.Infeasible)
		{
			numinf++;
			break;
		}

		
		//------------------get the optimal dual values-------
		 
		price[0]=MP.getDual(rng_3[0]);
		
		for (j=1;j<Jmax;j++)
		{ 
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}

		for (j=Jmax;j<depth;j++)
		{
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}
		 /*outfile << "Heuristic2 CG" << endl;
		 outfile << "Problem status: " << MP.getStatus() << endl;
		 outfile << "Objective of the RMP: " << MP.getObjValue() << endl;
		 outfile << "Values of decsision variables:" << endl;

		 for ( s = 0; s < decs.getSize(); s++) {
			 outfile << "lambda" << s << " = " << MP.getValue(decs[s]) <<" reduced cost = "<<MP.getReducedCost(decs[s])<< endl;
		 }
		 
		 outfile << "Prices:" <<endl;

		 for (j = 0; j < depth; j++) {
			  outfile << "con_" << j << " = " << price[j]<< endl; 
		 }*/

		 if(MP.getObjValue()<BestLB)
		 {  
			  imp_amount_h2=imp_amount_h2+(BestLB-MP.getObjValue());
			  h2imp++;
			  gear=1;
			  bayrak=1;
			  BestLB=MP.getObjValue(); 
			  //RMP_sim<<BestLB<<" 2"<<endl;
			  break;	
		 }
		else //break;
			Q++;

		if(Q>25) break;
		 

	}// end of heuristic pp2
	   
		if(bayrak==1)
			gear=1;
		else
			gear++;
	 }	//gear2 end
		
	 
     if(gear==3)
	 {
		bayrak=0;
		kol_say=0;
		for (;;)
		{

		reduced_cost.clear();
		reduced_cost-=price[0];
		for( j=1; j< Jmax; j++)
		{
			reduced_cost-=price[j]*y[j];
		}
		
		// adding new terms coming from newly added constraints
		
		for (j=Jmax;j<depth;j++)
		{
			if (hist[j-Jmax+1][0]==1)	
				reduced_cost-=price[j]*(y[hist[j-Jmax+1][1]]+y[hist[j-Jmax+1][2]]);
			else
				reduced_cost+=price[j]*(y[hist[j-Jmax+1][1]]+y[hist[j-Jmax+1][2]]);
		}

		
		// reduced problem pricing obj-------------------
		for( i=1; i< Imax; i++)
		 { 
			 for (j=1;j<=rn[i];j++)
		 		reduced_cost += x[i][order[i][j]]*d[i][order[i][j]] * h[i];		         				
		 }
		//--------------------------------------------
		obj1.setExpr(reduced_cost);

		numofcalle++;
		pptimer.restart();
		SP.solve();
		pptime=pptime+pptimer.stop();

		//off<<"pricing status: "<<SP.getStatus()<<endl;
		
		// -------------get the new column------------------------
		new_col[0]=1;
		
		for(j=1;j<Jmax;j++) 
		{
			if (SP.getValue(y[j]) < RC_EPS)
				new_col[j] = 0;
			else
				new_col[j] = 1;
		}
				
		for (j=Jmax;j<depth;j++) 
		{
			flag1=0;flag2=0;
			if (new_col[hist[j-Jmax+1][1]] ==1 )
				flag1=1;
			if (new_col[hist[j-Jmax+1][2]] ==1 )
				flag2=1;
		 
		    if(flag1==1 && flag2==1)
			   new_col[j]=1;
			else
			   new_col[j]=0;  
		   
		}
		cost=0;
		
		// for reduced pricing-----------------------
		for (i=1;i<Imax;i++)
		{
			min_dist=M;
			for (j=1;j<=rn[i];j++)
			{
				if (new_col[order[i][j]]==1 && d[i][order[i][j]]< min_dist)
				{
					min_distj=order[i][j];
					min_dist=d[i][order[i][j]];
				}
			}
			cost=cost+d[i][min_distj]*h[i];
		}
		//--------------------
		//---------------------------report_SP------------------------------
		
		/*outfile2 << "Problem status: " << SP.getStatus() << endl;
		outfile2 << "Objective of the PP: " << SP.getObjValue() << endl;
	    outfile2 << "New column:" << endl;
		for(j=1;j<Jmax;j++) 
		{
			outfile2 << "y" << j << " = " << SP.getValue(y[j]) << endl;
		}

		outfile2 << cost << endl;*/

		//off<<SP.getObjValue()<<endl;
		if (SP.getObjValue() > -RC_EPS) 
		{
			break;
		}

		numcg++;
		av_single_day_e=av_single_day_e+cost;
		//single_day<<cost<<" 3"<<endl;

		decs.add(IloNumVar(obj_exp(cost) + rng_3(new_col) ));
		kol_say++;
		index=kolon.getSize();
		kolon.add(IloNumArray (env, P+1)); 
		kolon[index][0]=cost;
		kolon[index][1]=1;
		for (j=1;j<Jmax;j++)
		    kolon[index][j+1]=new_col[j];

		mptimer.restart();
		MP.solve();
		mptime=mptime+mptimer.stop();
		//cpu_mp<<mptimer.stop()<<" "<<mptime<<endl;

		if (MP.getStatus()== MP.Infeasible)
		{
			numinf++;
			break;
		}
		
		//------------------get the optimal dual values-------
		 
		price[0]=MP.getDual(rng_3[0]);
		
		for (j=1;j<Jmax;j++)
		{ 
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}

		for (j=Jmax;j<depth;j++)
		{
			price[j]=MP.getDual(rng_3[j]);
			new_col[j]=0;
		}
		 
		 /*outfile << "Problem status: " << MP.getStatus() << endl;
		 outfile << "Objective of the RMP: " << MP.getObjValue() << endl;
		 outfile << "Values of decsision variables:" << endl;

		 for ( s = 0; s < decs.getSize(); s++) {
			outfile << "lambda" << s << " = " << MP.getValue(decs[s]) <<" reduced cost = "<<MP.getReducedCost(decs[s])<< endl;
		 }
		 outfile << "Prices:" <<endl;

		 for (j = 0; j < depth; j++) {
			  outfile << "con_" << j << " = " << price[j]<< endl; 
		 }*/

		if(MP.getObjValue()<BestLB)
		 {
			   imp_amount_e=imp_amount_e+(BestLB-MP.getObjValue());
			   exactimp++;
			   //gear=1;
			   //bayrak=1;
			   BestLB=MP.getObjValue();
			   //RMP_sim<<BestLB<<" 3"<<endl;
			   //break;
		 }
	}// end of exact pp
	  	// CG iteration
	  if(bayrak==1 ) 
	  		gear=1;
	  else
		  gear=4;
	
	}//gear3 end
//} //while loop end
} // infeasible end	
	col_debug<<kol_say<<endl;
	//---------------------------------------------------------end of the node------------------------------------------------------
	cputime=timer.stop();
	timer.start();
	
	if (MP.getStatus()!= MP.Infeasible && cputime < 9000)
	{
		say1++;
		branch_timer.restart();

		sch_timer.restart();
		sayac=0;
		for (s=0;s<kolon.getSize();s++)
		{
			if(MP.getValue(decs[s])> RC_EPS)
			{
			sch[0][sayac]=MP.getValue(decs[s]);
		
			for (j=1;j<Jmax;j++)
			{		
				if (kolon[s][j+1]==1) 
				sch[r[j]][sayac]= j;
			}
			sayac++;
			}
		}
		sch_time=sch_time+sch_timer.stop();

		
		isFractional=0;
		isBetter=0;

		
		// check if there exist some fractional schedules
		for (s=0;s<decs.getSize();s++)	
		{
			if (MP.getValue(decs[s]) > RC_EPS)
			{
				if (integer(MP.getValue(decs[s])) == 0)
					isFractional=1;
			}
		}	

		// compare optimal RMP of the parent node with the BestSch
		if(MP.getObjValue() - BestValue < -RC_EPS)
			isBetter=1;

		branch_time1=branch_time1 + branch_timer.stop();	

		if(isBetter==1)
		{
			if(isFractional==0)
			{
				
				say2++;
				branch_timer.restart();

				//-------------------------------------------------
				BestValue=MP.getObjValue();
				//UB_sim<<BestValue<<endl;
				t=0;
				for (s=0;s<sayac;s++)
				{
					while(sch[0][s]>1)
					{
						t++;
						for(k=1;k<Kmax;k++)
						{
							BestSch[t][k]=sch[k][s];
						}
						sch[0][s]=sch[0][s]-1;
					}
				}
				numcol=sayac;
				env2.end();
				//------------------------------------------------------------------
				branch_time2=branch_time2 + branch_timer.stop();	
			}
			else
			{
				say3++;
				branch_timer.restart();
				//-----------lets find a heuristic UB at this point-----------------
				UBvalue=UB(sch,sayac);
				if(UBvalue<BestValue)
				{
					UB_imp++;
					//off<<"BestSol is updated from "<<BestValue<<" to "<<UBvalue<<endl;
					BestValue=UBvalue;

					for (k=1;k<Kmax;k++)
					{
						for (s=1;s<Tmax;s++)	
							BestSch[k][s]=UB_sch[k][s];
					}
					//UB_sim<<BestValue<<endl;
				}
				
				//------------------------------------------------------------------
				
				// calculate together_jp values: how many times j and p are together
				for(j=1;j<Jmax;j++)
				{
					for(p=1;p<Jmax;p++)
						together[j][p]=0;
				}
				for (s=0;s<sayac;s++)
				{
					for(p=1;p<K;p++)
					{
						ecz1=sch[p][s];
						for (j=p+1;j<Kmax;j++)
						{
							ecz2=sch[j][s];
							if (ecz2<ecz1)
								together[ecz2][ecz1]=together[ecz2][ecz1]+sch[0][s];
							else
								together[ecz1][ecz2]=together[ecz1][ecz2]+sch[0][s];
						}
					}		
				}
				env2.end();
				// lets make the branching	
				IloNumArray brvar(env,3);
				min_dev=M;
			    for(j=1;j<Jmax;j++)
			    {
				   for(p=j+1;p<Jmax;p++)
				   {
					   if(integer(together[j][p]) == 0 && (abs(ceil(together[j][p])-together[j][p]-0.5) < min_dev))
					   {
						   min_dev=abs(ceil(together[j][p])-together[j][p]-0.5);
						   brvar[0]=j;
						   brvar[1]=p;
						   brvar[2]=together[j][p];
					  }
				   }
			    }
				
				branch_time3=branch_time3 + branch_timer.stop();	
				branch(hist,brvar);
				brvar.clear();brvar.end();
				
			}
			hist.clear();hist.end();

		}//isBetter if statement end 
		else
		{
			// here is the rigt point to free the memory-------------------------
			say4++;
			branch_timer.restart();
			env2.end();
			//------------------------------------------------------------------
			branch_time4=branch_time4 + branch_timer.stop();	
		}
		
	}//infeasibility if statement end
	else
	{
	say5++;
	branch_timer.restart();
	env2.end();
	branch_time5=branch_time5 + branch_timer.stop();

	}
	

}// end of for loop

}

void generate_initcols()
{
	
	int time,min_dist,min_distj;

	for (t=1;t<Tmax;t++)
	{
		for (j=1;j<Jmax;j++)
			init_y[j][t]=0;
	}

	for (k=1;k<Kmax;k++)
	{
		time=1;
		for (j=1;j<Jmax;j++)
		{
			if (r[j]==k)
			{
				for (t=time;t<time+n[j];t++)
				{
					init_y[j][t]=1;
				}
				time=time+n[j];				
			}
		}
	}
	/*for(j=1;j<Jmax;j++)
	{
		for(t=1;t<Tmax;t++)
			off<<init_y[j][t]<<" ";
		off<<endl;
	}*/
	
	//-------------find objective coefficients -----------------
	for (t=1;t<Tmax;t++)
	{
		init_obj[t]=0;
		for (i=1;i<Imax;i++)
		{
			min_dist=M;
			for (j=1;j<Jmax;j++)
			{
				if (init_y[j][t]==1 && d[i][j]< min_dist)
				{
					min_distj=j;
					min_dist=d[i][j];
				}
			}
			
			init_obj[t]=init_obj[t]+d[i][min_distj]*h[i];
		}
	}
}

void readInput(char s[])
{
    printf("readInput\n");
    
    FILE *finput;
	finput=fopen(s,"r");

	fscanf(finput,"%d",&I);
	fscanf(finput,"%d",&J);
	fscanf(finput,"%d",&T);
	fscanf(finput,"%d",&K);

	

    for (i=1 ;i<=I ;i++)
	{
		for (j=1 ;j<=J ;j++)
		{
			fscanf(finput,"%d",&d[i][j]);
		}
	}

    for (i=1 ;i<=I ;i++)
	{
		fscanf(finput,"%d",&h[i]);
    }

    for (j=1 ;j<=J ;j++)
	{
		fscanf(finput,"%d",&r[j]);
    }

    for (j=1 ;j<=J ;j++)
	{
		fscanf(finput,"%d",&n[j]);
    }
  	fclose(finput);
}

bool   integer(float k)
{                 
 int_timer.restart();
 if (k==0)  return true;  
 if (k>0)  return integer (k-1);
 int_time=int_time+int_timer.stop();
 return false;
}

void region_order()
{
   
	double app_dist[Jmax][Jmax];

    int i,j,jj,k,kk,Enkisa,Enorta;
    //int app_dist[Jmax][Jmax];
    int reg_dist[Kmax][Kmax];
    int reg_metric[Kmax];

    /*for(k=1;k<Kmax;k++)
        ordered_regs[k]=k;*/

    for(j=1;j<Jmax;j++)
    {
        for (jj=j+1;jj<Jmax;jj++)
        {
            Enkisa=M;
            for(i=1;i<Imax;i++)
            {
                if((d[i][j]+d[i][jj])<Enkisa)
                {
                    Enkisa=d[i][j]+d[i][jj];
                }
            }
            app_dist[j][jj]=Enkisa;
            app_dist[jj][j]=Enkisa;
        }
    }

    for(k=1;k<Kmax;k++)
    {
        for(kk=1;kk<Kmax;kk++)
            reg_dist[k][kk]=0;
    }

    for(k=1;k<Kmax;k++)
    {
        for(kk=k+1;kk<Kmax;kk++)
        {
            Enkisa=M;
            for(j=1;j<=pnumber[k];j++)
            {
                for(jj=1;jj<=pnumber[kk];jj++)
                {
                    if(app_dist[region_pset[k][j]][region_pset[kk][jj]]<Enkisa)
                        Enkisa=app_dist[region_pset[k][j]][region_pset[kk][jj]];
                }
            }
            reg_dist[k][kk]=Enkisa;
            reg_dist[kk][k]=Enkisa;
        }
    }

    //fprintf(fout,"bolge uzaklik matrisi\n");
    for(k=1;k<Kmax;k++)
    {
        reg_metric[k]=0;
        for(kk=1;kk<Kmax;kk++)
        {
            reg_metric[k]=reg_metric[k]+reg_dist[k][kk];
            //fprintf(fout,"%3d ",reg_dist[k][kk]);
        }

        //fprintf(fout,"\n");
    }

    for(k=1;k<Kmax;k++)
    {
        Enkisa=M;//0
        Enorta=0;
        for(kk=1;kk<Kmax;kk++)
        {
            //fprintf(fout,"%d ",reg_metric[kk]);
            if(reg_metric[kk]<Enkisa)//>
            {
                Enkisa=reg_metric[kk];
                Enorta=kk;
            }
        }
        ordered_regs[k]=Enorta;
        reg_metric[Enorta]=M;//0
        //fprintf(fout,"%d\n",ordered_regs[k]);
    }
    /*for(k=1;k<Kmax;k++)
        fprintf(fout,"%d ",ordered_regs[k]);
    fprintf(fout,"\n");*/
    
}

int heuristic1 (IloNumArray fiyat, IloNumArray ciz)
{

	int i,j,k;
	int min_k;
	double mindist_k;
	int maliyet;

	

	for (j=1;j<Jmax;j++)
	{
		ciz[j]=0;
	}

	for(j=1;j<Jmax;j++)
	{
		if(pd[j]!=0)
			avg[j]=(double)(phar_dset[j][pd[j]+1]-fiyat[j])/phar_dset[j][pd[j]+2];
		else
			avg[j]=M;
	}

	// writing the entries
	/*off<<"exact assignment costs"<<endl;
	for (j=1;j<Jmax;j++)
	{
		for (i=1;i<=pd[j]+2;i++)
			off<<phar_dset[j][i]<<" ";
		off<<avg[j]<<endl;
	}*/

	for (k=1;k<Kmax;k++)
	{
		mindist_k=M;
		for(j=1;j<=pnumber[k];j++)
		{
			if(avg[region_pset[k][j]]<mindist_k)
			{
				mindist_k=avg[region_pset[k][j]];
				min_k=region_pset[k][j];
			}
		}
		if(mindist_k==M)
			ciz[region_pset[k][1]]=1;
		else
			ciz[min_k]=1;
	}

	// calculate the exact assignment cost
	maliyet=0;

	for (i=1;i<Imax;i++)
	{
		min_dist=M;
		for (j=1;j<Jmax;j++)
		{
			if (ciz[j]==1 && d[i][j]< min_dist)
			{
				min_distj=j;
				min_dist=d[i][j];
			}
		}
		maliyet=maliyet+d[i][min_distj]*h[i];
	}
	
	return(maliyet);
}

void find_exactA()
{
	int i,j;

	for(j=1;j<Jmax;j++)
		pd[j]=0;

	for(i=1;i<Imax;i++)
	{
		pd[order[i][1]]=pd[order[i][1]]+1;
		phar_dset[order[i][1]][pd[order[i][1]]]=i;
		

		for(j=2;j<=rn[i];j++)
		{
			if(r[order[i][j]] == r[order[i][1]])
			{
				pd[order[i][j]]=pd[order[i][j]]+1;
				phar_dset[order[i][j]][pd[order[i][j]]]=i;
			}
			else
				break;
		
		}
	
	}

	for (j=1;j<Jmax;j++)
	{
		phar_dset[j][pd[j]+1]=0;
		phar_dset[j][pd[j]+2]=0;
		for (i=1;i<=pd[j];i++)
		{
			phar_dset[j][pd[j]+1]=phar_dset[j][pd[j]+1]+h[phar_dset[j][i]]*d[phar_dset[j][i]][j];
			phar_dset[j][pd[j]+2]=phar_dset[j][pd[j]+2]+h[phar_dset[j][i]];
		}
		
	}

}

int heuristic2 (IloNumArray fiyat, IloNumArray ciz)
{
	int i,j,k,current_k,current_p,min_dist,min_distj,j2;
	double current_cost,best_cost,best_p;
	int maliyet;

	for (j=1;j<Jmax;j++)
	{
		ciz[j]=0;
	}

	//off<<"heuristic3"<<endl;

	
	for (k=1;k<Kmax;k++)
	{
		best_cost=M;
		current_k=ordered_regs[k];
		
		for (j=1;j<=pnumber[current_k];j++)
		{
			current_cost=0;
			current_p=region_pset[current_k][j];

			// calcualte the assignment cost
			ciz[current_p]=1;
			maliyet=0;
			for (i=1;i<Imax;i++)
			{
				min_dist=M;
				for (j2=1;j2<Jmax;j2++)
				{
					if (ciz[j2]==1 && d[i][j2]< min_dist)
					{
						min_distj=j2;
						min_dist=d[i][j2];
					}
				}		
				maliyet=maliyet+d[i][min_distj]*h[i];
			}

			current_cost=maliyet;
			for (j2=1;j2<Jmax;j2++)
			{
				if(ciz[j2]==1)
					current_cost=current_cost-fiyat[j2];
			}

			ciz[current_p]=0;
			
			if(current_cost<best_cost)
			{
				best_cost=current_cost;
				best_p=current_p;
			}
			//off<<current_k<<" "<<current_p<<" "<<current_cost<<endl;
		}
		ciz[best_p]=1;
		//off<<best_cost<<" "<<best_p<<endl;
	}
	

	// calculate the exact assignment cost
	maliyet=0;
	for (i=1;i<Imax;i++)
	{
		min_dist=M;
		for (j=1;j<Jmax;j++)
		{
			if (ciz[j]==1 && d[i][j]< min_dist)
			{
				min_distj=j;
				min_dist=d[i][j];
			}
		}
		maliyet=maliyet+d[i][min_distj]*h[i];
	}

	return(maliyet);
}

int UB(IloNumArrayArray frac_sch, int boyut )
{ 

	IloEnv env3;
	int i,j,k,t,s,Tint,Trem,EnYakin,EnYakinj;
	int temp_n[Jmax];

	IloNumArrayArray temp_sch(env3, Kmax);
	for(k=0; k < Kmax; k++)
    {
		temp_sch[k] = IloNumArray (env3, sayac);   
    }
	
				
	for(k=0;k<Kmax;k++)
	{				
		for (s=0;s<boyut;s++)
			temp_sch[k][s]=frac_sch[k][s];
	}
	

	for(j=1;j<Jmax;j++)
		temp_n[j]=n[j];

	int maliyet;

	/*off<<"frac_sch:"<<endl;
	// frac_sch yi yazalim------------------------------
	for(k=0;k<Kmax;k++)
	{				
		for (s=0;s<boyut;s++)
			off<<frac_sch[k][s]<<"   ";
		off<<endl;
	}*/
	//---------------------------------------------

	t=0;
	for (s=0;s<boyut;s++)
	{
		while(temp_sch[0][s]>=1)
		{
			t++;
			for(k=1;k<Kmax;k++)
				UB_sch[k][t]=temp_sch[k][s];
			temp_sch[0][s]=temp_sch[0][s]-1;	
		}
	}
	Tint=t;

	
	ToTint=ToTint+Tint;
	CoTint++;
	Trem=Tmax-1-Tint;

	if (Tint >= 0) //(Tmax-1)*0.75
	{
	noUsed++;
	// update the duties and run a random UB algorithm
	for(t=1;t<=Tint;t++)
	{
		for(k=1;k<Kmax;k++)
		{
			temp_n[UB_sch[k][t]]=temp_n[UB_sch[k][t]]-1;
		}
	}
	for(t=Tint+1;t<Tmax;t++)
	{
		for(k=1;k<Kmax;k++)
		{
			for(j=1;j<=pnumber[k];j++)
			{
				if(temp_n[region_pset[k][j]]>0)
				{
					UB_sch[k][t]=region_pset[k][j];
					temp_n[region_pset[k][j]]=temp_n[region_pset[k][j]]-1;
					j=pnumber[k]+1;
				}
			}
		}
	}

maliyet=0;

for (i=1;i<Imax;i++)
    {
         for (t=1;t<Tmax;t++)
         {
              EnYakin=M;
              for(k=1;k<Kmax;k++)
              {
                  if(d[i][UB_sch[k][t]]< EnYakin)
                  {
                  EnYakin= d[i][UB_sch[k][t]];
                  EnYakinj=UB_sch[k][t];
                  }
              }
              maliyet=maliyet+EnYakin*h[i];
         }
     }
	/*off<<"UB_sch with cost: "<<maliyet<<endl;
	// UBsch yi yazalim------------------------------
	for(k=1;k<Kmax;k++)
	{				
		for (s=1;s<Tmax;s++)
			off<<UB_sch[k][s]<<"   ";
		off<<endl;
	}*/
	//---------------------------------------------
	}
	else
		maliyet=M;
env3.end();
return(maliyet);


}

int node_heuristic1 (IloNumArray fiyat, IloNumArrayArray ek, IloNumArray ciz)
{

	int i,j,k,l;
	int min_k;
	double mindist_k;
	int maliyet;

	for (j=1;j<Jmax;j++)
	{
		ciz[j]=0;
	}

	for(j=1;j<Jmax;j++)
	{
		for (l=1;l<ek.getSize();l++)
		{
			if(ek[l][1]==j || ek[l][2]==j)
			{
				if(ek[l][0]==1)
					phar_dset[j][pd[j]+1]=phar_dset[j][pd[j]+1]-fiyat[Jmax+l-1];
				else
					phar_dset[j][pd[j]+1]=phar_dset[j][pd[j]+1]+fiyat[Jmax+l-1];
			}	  
		}
	
		if(pd[j]!=0)
			avg[j]=(double)(phar_dset[j][pd[j]+1]-fiyat[j])/phar_dset[j][pd[j]+2];
		else
			avg[j]=M;
	}

	// writing the entries
	/*off<<"exact assignment costs"<<endl;
	for (j=1;j<Jmax;j++)
	{
		for (i=1;i<=pd[j]+2;i++)
			off<<phar_dset[j][i]<<" ";
		off<<avg[j]<<endl;
	}*/

	for (k=1;k<Kmax;k++)
	{
		mindist_k=M;
		for(j=1;j<=pnumber[k];j++)
		{
			if(avg[region_pset[k][j]]<mindist_k)
			{
				mindist_k=avg[region_pset[k][j]];
				min_k=region_pset[k][j];
			}
		}
		if(mindist_k==M)
			ciz[region_pset[k][1]]=1;
		else
			ciz[min_k]=1;
	}

	// calculate the exact assignment cost
	maliyet=0;

	for (i=1;i<Imax;i++)
	{
		min_dist=M;
		for (j=1;j<Jmax;j++)
		{
			if (ciz[j]==1 && d[i][j]< min_dist)
			{
				min_distj=j;
				min_dist=d[i][j];
			}
		}
		maliyet=maliyet+d[i][min_distj]*h[i];
	}
	
	return(maliyet);
	
			
}

int node_heuristic2 (IloNumArray fiyat, IloNumArrayArray ek, IloNumArray ciz)
{
	int i,j,k,current_k,current_p,min_dist,min_distj,j2;
	double current_cost,best_cost,best_p;
	int maliyet;

	for (j=1;j<Jmax;j++)
	{
		ciz[j]=0;
	}

	//off<<"heuristic3"<<endl;

	
	for (k=1;k<Kmax;k++)
	{
		best_cost=M;
		current_k=ordered_regs[k];
		
		for (j=1;j<=pnumber[current_k];j++)
		{
			current_cost=0;
			current_p=region_pset[current_k][j];

			// calcualte the assignment cost
			ciz[current_p]=1;
			maliyet=0;
			for (i=1;i<Imax;i++)
			{
				min_dist=M;
				for (j2=1;j2<Jmax;j2++)
				{
					if (ciz[j2]==1 && d[i][j2]< min_dist)
					{
						min_distj=j2;
						min_dist=d[i][j2];
					}
				}		
				maliyet=maliyet+d[i][min_distj]*h[i];
			}

			current_cost=maliyet;
			for (j2=1;j2<Jmax;j2++)
			{
				if(ciz[j2]==1)
				{
					current_cost=current_cost-fiyat[j2];
					for (l=1;l<ek.getSize();l++)
					{
						if(ek[l][1]==ciz[j2] || ek[l][2]==ciz[j2])
						{
							if(ek[l][0]==1)
								current_cost=current_cost-fiyat[Jmax+l-1];
							else
								current_cost=current_cost+fiyat[Jmax+l-1];
						}	  
					}
				}
				
			}

			ciz[current_p]=0;
			
			if(current_cost<best_cost)
			{
				best_cost=current_cost;
				best_p=current_p;
			}
			//off<<current_k<<" "<<current_p<<" "<<current_cost<<endl;
		}
		ciz[best_p]=1;
		//off<<best_cost<<" "<<best_p<<endl;
	}
	

	// calculate the exact assignment cost
	maliyet=0;
	for (i=1;i<Imax;i++)
	{
		min_dist=M;
		for (j=1;j<Jmax;j++)
		{
			if (ciz[j]==1 && d[i][j]< min_dist)
			{
				min_distj=j;
				min_dist=d[i][j];
			}
		}
		maliyet=maliyet+d[i][min_distj]*h[i];
	}

	return(maliyet);
}

void heuristic1_analysis()
{
	// we use 6090159_6 instance for this analysis.
	int i,j,k,m;
	int min_k;
	double mindist_k;
	int maliyet;
	double fiyat[Jmax];
	int ciz[Jmax];

	FILE *price_input,*h1_analysis;
	price_input=fopen("prices.txt","r");
	h1_analysis=fopen("pp_h1_results.txt","w");

	

	for (m=1;m<=573;m++)
	{
		for(j=0;j<Jmax;j++)
		{
			fscanf(price_input,"%lf",&fiyat[j]);
			ciz[j]=0;
		}
	
		h1_timer.restart();
		for(j=1;j<Jmax;j++)
		{
			if(pd[j]!=0)
				avg[j]=(double)(phar_dset[j][pd[j]+1]-fiyat[j])/phar_dset[j][pd[j]+2];
			else
				avg[j]=M;
		}

		for (k=1;k<Kmax;k++)
		{
			mindist_k=M;
			for(j=1;j<=pnumber[k];j++)
			{
				if(avg[region_pset[k][j]]<mindist_k)
				{
					mindist_k=avg[region_pset[k][j]];
					min_k=region_pset[k][j];
				}
			}
			if(mindist_k==M)
				ciz[region_pset[k][1]]=1;
			else
				ciz[min_k]=1;
		}

		// calculate the exact assignment cost
		maliyet=0;

		for (i=1;i<Imax;i++)
		{
			min_dist=M;
			for (j=1;j<Jmax;j++)
			{
				if (ciz[j]==1 && d[i][j]< min_dist)
				{
					min_distj=j;
					min_dist=d[i][j];
				}
			}
			maliyet=maliyet+d[i][min_distj]*h[i];
		}

		maliyet=maliyet-fiyat[0];
		for (j=1;j<Jmax;j++) maliyet=maliyet-fiyat[j]*ciz[j];

		fprintf(h1_analysis,"%.10f %d",h1_timer.stop(),maliyet);
		fprintf(h1_analysis,"\n");
		//h1_analysis << h1_timer.stop() << " " << maliyet << endl;

	}

	fclose(h1_analysis);

}

void heuristic2_analysis()
{
	// we use 6090159_6 instance for this analysis.
	int m;
	int i,j,k,current_k,current_p,min_dist,min_distj,j2,best_p;
	double current_cost,best_cost;
	int maliyet;
	double fiyat[Jmax];
	int ciz[Jmax];

	FILE *price_input,*h2_analysis;
	price_input=fopen("prices.txt","r");
	h2_analysis=fopen("pp_h2_results.txt","w");

	

	for (m=1;m<=573;m++)
	{
		for(j=0;j<Jmax;j++)
		{
			fscanf(price_input,"%lf",&fiyat[j]);
			ciz[j]=0;
		}
	
		h2_timer.restart();

		for (k=1;k<Kmax;k++)
		{
			best_cost=M;
			current_k=ordered_regs[k];
		
			for (j=1;j<=pnumber[current_k];j++)
			{
				current_cost=0;
				current_p=region_pset[current_k][j];

				// calcualte the assignment cost
				ciz[current_p]=1;
				maliyet=0;
				for (i=1;i<Imax;i++)
				{
					min_dist=M;
					for (j2=1;j2<Jmax;j2++)
					{
						if (ciz[j2]==1 && d[i][j2]< min_dist)
						{
							min_distj=j2;
							min_dist=d[i][j2];
						}
					}		
					maliyet=maliyet+d[i][min_distj]*h[i];
				}

				current_cost=maliyet;
				for (j2=1;j2<Jmax;j2++)
				{
					if(ciz[j2]==1)
						current_cost=current_cost-fiyat[j2];
				}

				ciz[current_p]=0;
			
				if(current_cost<best_cost)
				{
					best_cost=current_cost;
					best_p=current_p;
				}
				//off<<current_k<<" "<<current_p<<" "<<current_cost<<endl;
			}
			ciz[best_p]=1;
			//off<<best_cost<<" "<<best_p<<endl;
		}
	

		// calculate the exact assignment cost
		maliyet=0;
		for (i=1;i<Imax;i++)
		{
			min_dist=M;
			for (j=1;j<Jmax;j++)
			{
				if (ciz[j]==1 && d[i][j]< min_dist)
				{
					min_distj=j;
					min_dist=d[i][j];
				}
			}
			maliyet=maliyet+d[i][min_distj]*h[i];
		}

		maliyet=maliyet-fiyat[0];
		for (j=1;j<Jmax;j++) maliyet=maliyet-fiyat[j]*ciz[j];

		fprintf(h2_analysis,"%.5f %d",h2_timer.stop(),maliyet);
		fprintf(h2_analysis,"\n");
		

	}

	fclose(h2_analysis);

}