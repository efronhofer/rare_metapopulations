{
									
	Simulations of Insects in Structured Populations			
			-----------------------------			
		 Field Station Fabrikschleichach		
			Evolutionary Ecology Group			
			  University of Wuerzburg			
									
	    	 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.
        
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
        
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
    MA 02110-1301, USA.
}



program SISP;

//______________________________________________________________________________
//---------------------------------------------------------------------Constants

const
  XMAX = 10;		//max amount of patches in x-direction
  YMAX = 10;		//max amount of patches in y-direction
  MAXF = 30000; 		//max amount of females per patch
  MAXM = 30000;		//max amount of males per patch
  VARIANCE = 0.15;	//variance used for calculation of many things
  MAXMUT = 0.05;	//maximum possible value, that is added by mutation to
			//to the dispersal propensity alleles
  NOALLELEVALUES = 100;

//______________________________________________________________________________
//-------------------------------------------------------------------------Types

type

TAnalysis = object

  occupied : Array[1..XMAX,1..YMAX] of boolean;		//is patch occupied? 1/0
  turnover : double;		    //turnover per x-stripe
  pop_b4_rep : Array[1..XMAX,1..YMAX] of integer;
  pop_after_rep : Array[1..XMAX,1..YMAX] of integer;	//pops. after reproduction
  pop_last_gen : Array[1..XMAX,1..YMAX] of integer;		//pops. at t-1 for sink calc.
  pop_last_gen_m : Array[1..XMAX,1..YMAX] of integer;	
  pop_last_gen_f : Array[1..XMAX,1..YMAX] of integer;	
  mean_disp : Array[1..XMAX,1..YMAX] of double;		//mean dispersal rate
  mean_disp_na : Array[1..XMAX,1..YMAX] of boolean;	//was this patch occupied
							//before dispersal?
  mean_disp_all : double;				//total mean disp rate
  meta_occupied : integer;				//total amount of occ.
  meta_occupied_pre : integer;
  meta_occupied_post : integer;
  meta_occ_pre : double;
  meta_occ_post : double;
  meta_occ : real;
  Fst : real;
  emigrants_all : int64;				//total nr of emis
  metapop_pre : int64;
  metapop_post : int64;

  colonized : integer;
  extinct : integer;
  sum_ind_b4_all : int64;				//total nr of inds b4
  count_loc1_all : int64;
  sum_loc1_all : double;
  
  rescued : integer;	//number of rescued patches
  resc : double;	//fraction of rescued patches
  
  sum_ind_b4 : Array[1..XMAX,1..YMAX] of integer;
  emigrants : Array[1..XMAX,1..YMAX] of integer;
  immigrants : Array[1..XMAX,1..YMAX] of integer;
  emigrants_m : Array[1..XMAX,1..YMAX] of integer;
  immigrants_m : Array[1..XMAX,1..YMAX] of integer;
  emigrants_f : Array[1..XMAX,1..YMAX] of integer;
  immigrants_f : Array[1..XMAX,1..YMAX] of integer;
  
  occ_after_emi : double;
  occ_after_immi : double;
  
  changes_pre : integer;
  changes_post : integer;
  turnover_pre : double;
  turnover_post : double;
  rescueds : integer;

end;

TInd = object

  loc1_allel1,loc1_allel2 : real;	//2 alleles at locus 1 for
  					//calculation of dispersal trait
  loc2 : Array[1..2] of integer;
					//neutral loci for Fst-calculation
  dispersed : boolean;			//already dispersed or not?
  procedure InitInd;			//initialize individual randomly
  procedure genetics (var kid: TInd);
  						//creates offspring with random
  						//sex depending on the parents
  						//alleles
end;

TPatch = object

  fpopsize : integer;			//female population size
  mpopsize : integer;			//male population size
  popsize : integer;			//total population size
  male : Array[1..MAXM] of TInd;	//males of population
  female : Array[1..MAXF] of TInd;	//females of population
  capacity : integer;			//habitat capacity
  sigma : real;				//environmental fluctuations
  lambda_0 : real;			//per capita growth rate per individual
  disp_mort : real;			//dispersal mortality
  ext_prob : real;			//extinction probability
  occ_post : boolean;
  occ_pre : boolean;
  occ_pre_b4 : boolean;
  occ_post_b4 : boolean;

end;

TWorld = Array[1..XMAX,1..YMAX] of TPatch;	//whole metapopulation

TParameters = object

  capacity : integer;		//habitat capacity to be read in
  sigma : real;			//environmental fluctuations to be read in
  lambda_0 : real;		//per capita growth rate to be read in
  ext_prob : real;		//extinction probability
  disp_mort : real;		//dispersal mortality to be read in
  disp_start : real;		//starting value of traits
  ddd : boolean;		//density-dependent or independent disp.?
  allee : real;			//critical density for allee effect
  mut_prob : real;		//mutation probability of genetic traits
  tmax : integer;		//maximum no. of generations to be simulated
  max_x,
  max_y : integer;		//world dimensions
  grad_check : Array[1..5] of boolean;
  				//saves, what gradient should be chosen:
  				//1 - k
  				//2 - lambda_0
  				//3 - sigma
  				//4 - disp_mort
  				//5 - ext_prob

  max_capacity : integer;       //maximum capacity occuring in the world

  y_1,y_2 : real;		//first and second values of gradient
  slope : real;                 //slope of the sigmoid gradient
  position : integer;           //position (added / subtracted from worlds
                                //middle) of the inflection point of the grad

  gradient : boolean;		//if a gradient is used

  sigmoid : boolean;            //sigmoid gradient?

  range_expansion : boolean;	//are individuals brought out at x=x/2 or in the
  				//in the whole world?

  replications : integer;	//how often is one experiment replicated?

  offspr_exceeded : integer;	//how often was the maximum offspring number
  				//exceeded?

  procedure set_parameters(var infile : textfile);
  				//setting all parameters by an infile

end;

//______________________________________________________________________________
//---------------------------------------------------------------------Variables

var
  world : TWorld;		//meta-population
  ana : TAnalysis;		//all values necessary for analysis
  par : TParameters;		//all simulation parameters
  xsource,
  ysource,
  xtarget,
  ytarget : integer;		//running variables for dispersal matrix

  density : Array[1..XMAX,1..YMAX] of double;
					//patch densities for density-dependent dispersal

  outfile : textfile;

  parinfile : textfile;
 

//______________________________________________________________________________
//----------------------------------------------------------------Initialization

procedure init;
begin
  
  assign(outfile,'data/output.txt');
  rewrite(outfile);
  writeln(outfile,
't	K   sigma  mort   lamb  rep	Fst      Ipre	Ipost        Tpre	Tpost       emi	R');
  assign(parinfile,'data/para.in');
  reset(parinfile);

  par.set_parameters(parinfile);
end;

//______________________________________________________________________________
//--------------------------------------------------------Distribution Functions

function Poisson(lamb:double):integer;
// this function creates Poisson distributed random numbers
// with mean lamb. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
end;

function gauss : real;
// this function creates Gaussian distributed random numbers with mean
// 0 and standard deviation 1. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
 end;

function loggauss(Fertility, sigma: real): real;
// calculates the mean random fertility
// of a single female individual from a
// Gaussian probability distribution,
// using the variable "Fertility" as
// mean value and "sigma" (environmental
// fluctuations) as standard deviation. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
end;

//______________________________________________________________________________
//----------------------------------------------------------------------Analysis

procedure clear_analysis_disp;	// reset all analysis structures
var
  x,y : integer;
begin

  ana.emigrants_all:=0;
  ana.sum_ind_b4_all:=0;
  ana.sum_loc1_all:=0;
  ana.count_loc1_all:=0;
  ana.extinct:=0;
  ana.colonized:=0;
  
  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
  begin
  ana.emigrants[x,y]:=0;
  ana.immigrants[x,y]:=0;
  ana.emigrants_m[x,y]:=0;
  ana.immigrants_m[x,y]:=0;
  ana.emigrants_f[x,y]:=0;
  ana.immigrants_f[x,y]:=0;
  end;

end;

procedure analyze_disp;	// calculate emigration rates
begin

  if ana.sum_ind_b4_all>0 then ana.mean_disp_all:=ana.emigrants_all/ana.sum_ind_b4_all
  else ana.mean_disp_all:=0;

end;

function calculate_Fst : double;	// calculate FST-values
var
  x,y,f,m,i : integer;
  number : int64;
  p_hat : Array[1..XMAX,1..YMAX,1..NOALLELEVALUES] of double;
  p_hat_sum : Array[1..NOALLELEVALUES] of double;
  count : integer;
  Hs_sum_patch : Array[1..XMAX,1..YMAX] of double;
  Hs : double;
  Ht : double;
  Hs_sum : double;
  p_bar : Array[1..NOALLELEVALUES] of double;
  Ht_sum : double;
begin

  for i:=1 to NOALLELEVALUES do p_hat_sum[i]:=0;
  
  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
  begin
  Hs_sum_patch[x,y]:=0;
  end;
  
  count:=0;
  Hs_sum := 0;
  Ht_sum := 0;

  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
  begin
  
    if (world[x,y].fpopsize+world[x,y].mpopsize)>0 then
    begin

    for i:=1 to NOALLELEVALUES do
    begin
  
      number:=0;
  
      for f:=1 to world[x,y].fpopsize do
      begin
        if world[x,y].female[f].loc2[1]=i then inc(number);
        if world[x,y].female[f].loc2[2]=i then inc(number);
      end;
      for m:=1 to world[x,y].mpopsize do
      begin
        if world[x,y].male[m].loc2[1]=i then inc(number);
        if world[x,y].male[m].loc2[2]=i then inc(number);
      end;
    
      p_hat[x,y,i] := number/(2*(world[x,y].fpopsize+world[x,y].mpopsize));
      p_hat_sum[i]:=p_hat_sum[i]+p_hat[x,y,i];  
      Hs_sum_patch[x,y]:=Hs_sum_patch[x,y]+(p_hat[x,y,i]*p_hat[x,y,i]);
    end;
    inc(count);
    
    Hs_sum_patch[x,y]:=1-Hs_sum_patch[x,y];
    Hs_sum := Hs_sum+Hs_sum_patch[x,y];

    end;
  end;

  if count>0 then
  begin

  for i:=1 to NOALLELEVALUES do
  begin
    p_bar[i]:=p_hat_sum[i]/count;
    Ht_sum := Ht_sum+(p_bar[i]*p_bar[i]);   
  end;
  
  Hs:=Hs_sum/count;
  Ht:=1-Ht_sum;
  
  if Ht>0 then calculate_Fst := (Ht-Hs)/Ht else calculate_Fst:=0;
  
  end else calculate_Fst:=-3;

end;
  

procedure analyze (var outfile : textfile; const t,repli : integer);
						//calculation of mean
var 					//population density, dispersal rate of
  x,y : integer;		//occupied patches, portion of occupied
  c : double;			//index, that indicates source (c>1) or
						//sink (c<1) populations
  
  whole_sinks : double;
 
  sum_whole_sinks : integer;
  sosi_whole_count : integer;
   
begin

  analyze_disp;

  sum_whole_sinks:=0;
  sosi_whole_count:=0;
  ana.meta_occupied:=0;

  for x:=1 to par.max_x do
    begin //1
 
    for y:=1 to par.max_y do
      begin  //2


      if (world[x,y].mpopsize+world[x,y].fpopsize)>0 then
         inc(ana.meta_occupied);

      if (ana.pop_last_gen[x,y]>0) then
          begin
	  c:=(ana.emigrants[x,y]+ana.pop_after_rep[x,y]-ana.immigrants[x,y])/ana.pop_last_gen[x,y];
	  if (c<1) then
            begin
            inc(sum_whole_sinks);
            end;
      inc(sosi_whole_count);
	  end;

      if ((ana.pop_last_gen[x,y]=0) and (ana.immigrants[x,y]>0)) then
        begin
        inc(sum_whole_sinks);
        inc(sosi_whole_count);
        end;
        
      if ((((ana.pop_last_gen_m[x,y]-ana.emigrants_m[x,y])<1) or
          ((ana.pop_last_gen_f[x,y]-ana.emigrants_f[x,y])<1)) and
          (((ana.pop_last_gen_m[x,y]-ana.emigrants_m[x,y]+ana.immigrants_m[x,y])>=1) and
          ((ana.pop_last_gen_f[x,y]-ana.emigrants_f[x,y]+ana.immigrants_f[x,y])>=1)))then
        inc(ana.rescued);  

      end;  //2

    end; //1

    ana.meta_occ:=ana.meta_occupied/(par.max_x*par.max_y);
    ana.meta_occ_pre:=ana.meta_occupied_pre/(par.max_x*par.max_y);
    ana.meta_occ_post:=ana.meta_occupied_post/(par.max_x*par.max_y);

    if sosi_whole_count>0 then
      begin
      whole_sinks:=sum_whole_sinks/(sosi_whole_count)
      end else
      begin
      whole_sinks:=0;
      end;    

    //TURNOVER CALCULATIONS:
    
    ana.changes_pre:=0;
    ana.changes_post:=0;
    ana.rescueds:=0;
    
    for x:=1 to par.max_x do
    for y:=1 to par.max_y do
    begin
		if (world[x,y].occ_pre<>world[x,y].occ_pre_b4) then inc(ana.changes_pre);
		if (world[x,y].occ_post<>world[x,y].occ_post_b4) then inc(ana.changes_post);
		if ((world[x,y].occ_pre=false) and (world[x,y].occ_post=true)) then inc(ana.rescueds);
    end;
    
    
    ana.turnover_pre:=ana.changes_pre/(par.max_x*par.max_y);
    ana.turnover_post:=ana.changes_post/(par.max_x*par.max_y);
    ana.resc:=ana.rescueds/(par.max_x*par.max_y);
 
    ana.Fst:=calculate_Fst;

    writeln(outfile,t,'  ',par.capacity,'  ',par.sigma:5:3,'  ',par.disp_mort:5:3,
            '  ',par.lambda_0:5:3,'  ',repli,'  ',ana.Fst:7:5,'  ',
            ana.meta_occ_pre:7:5,'  ',ana.meta_occ_post:7:5,'  ', ana.turnover_pre:7:5,'  ',ana.turnover_post:7:5,'  ',
            '  ',ana.mean_disp_all:7:5,'  ',ana.resc:7:5)
	
end;			

//______________________________________________________________________________
//------------------------------------Calculate Occupancies b4/after immigration

procedure calc_occs;
var
  x,y,m,f : integer;
  mdisp,fdisp : integer;
  count_occ_after_emi,count_occ_after_immi : integer;
begin

  count_occ_after_emi:=0;
  count_occ_after_immi:=0;

  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
  begin
    mdisp:=0;
    fdisp:=0;       
    for m:=1 to world[x,y].mpopsize do
      if world[x,y].male[m].dispersed then inc(mdisp);
    for f:=1 to world[x,y].fpopsize do
      if world[x,y].female[f].dispersed then inc(fdisp);
        
    if (((world[x,y].fpopsize+world[x,y].mpopsize)-mdisp-fdisp)>0) then 
      inc(count_occ_after_emi);
    if (world[x,y].fpopsize+world[x,y].mpopsize)>0 then
      inc(count_occ_after_immi); 
  end;
  
  ana.occ_after_emi:=count_occ_after_emi/(par.max_x*par.max_y);
  ana.occ_after_immi:=count_occ_after_immi/(par.max_x*par.max_y);
  
end;
			
//______________________________________________________________________________
//------------------------------------------------------------Setting Parameters

procedure TParameters.set_parameters(var infile : textfile);
var
  d : string;
  rand : integer;		//random seed
begin

  readln(infile);
  readln(infile,d);
  readln(infile);
  readln(infile);
  readln(infile,max_x,max_y);
  readln(infile);
  readln(infile);
  readln(infile,capacity,sigma,lambda_0,disp_mort,
  		mut_prob,ext_prob,tmax,disp_start,allee,rand,
  		replications);
		
  RandSeed:=rand;
  
    if d='yes' 	then ddd:=true
				else ddd:=false;

  if ddd then disp_start:=0.9
		 else disp_start:=0.1;

end;

//______________________________________________________________________________
//---------------------------------------------------------Initialize Individual

procedure TInd.InitInd;
var
  a : integer;
begin

  loc1_allel1:=par.disp_start+VARIANCE*(random-0.5);
  loc1_allel2:=par.disp_start+VARIANCE*(random-0.5);
  
  for a:=1 to 2 do
  begin
    loc2[a]:=random(NOALLELEVALUES)+1;
  end;
  
  dispersed:=false;

end;

//______________________________________________________________________________
//-----------------------------------------------------Initialize Metapopulation

procedure makepop;
var
  x,y : integer;
  ind : integer;
  sex : real;
begin
  ana.meta_occupied:=0;
  par.offspr_exceeded:=0;

  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
    begin  //1

    with world[x,y] do
      begin  //2
      fpopsize:=0;
      mpopsize:=0;
      popsize:=0;
      capacity:=par.capacity;
      sigma:=par.sigma;
      lambda_0:=par.lambda_0;
      disp_mort:=par.disp_mort;
      ext_prob:=par.ext_prob;

        for ind:=1 to capacity do		//initialize every individual
          begin //4
          inc(popsize);
          inc(fpopsize);
          female[fpopsize].InitInd;		//female gets alleles
          end; //4

      end;  //2

    end;  //1

end;

//______________________________________________________________________________
//----------------------------------------------------------------------Genetics

 function mutation(const allele : real): real;
 begin
	mutation:=1+VARIANCE*gauss;
 end;

procedure TInd.genetics(var kid: TInd);
begin

  kid.dispersed:=false;
  kid.loc1_allel1:=loc1_allel1;
  kid.loc1_allel2:=loc1_allel2;
  kid.loc2[1]:=loc2[1];
  kid.loc2[2]:=loc2[2];

  if random<par.mut_prob then
                kid.loc1_allel1:=kid.loc1_allel1*
										mutation(kid.loc1_allel1);
  if random<par.mut_prob then
				kid.loc1_allel2:=kid.loc1_allel2*
										mutation(kid.loc1_allel2);

  if random<par.mut_prob then kid.loc2[1]:=random(NOALLELEVALUES)+1;
  if random<par.mut_prob then kid.loc2[2]:=random(NOALLELEVALUES)+1;


end;

//______________________________________________________________________________
//------------------------------------------------------------------Reproduction

// as this is the asexual version, only females are considered for
// reproduction and birth in the following. to simulate sexually
// reproducing organisms, mating and random birth of males needs to be
// added.

procedure reproduction(x,y : integer; lambda : real; 
		       newfpopsize,newmpopsize : integer);
var
  sqr_dens : double;			//squared pop density per patch
  a,
  survival : double;			//parameters for logistic growth
  sex : real;				//determines the sex of the newborns
  newfemale : Array[1..MAXF] of TInd;
  					//newborn female individuals
  newmale : Array[1..MAXM] of TInd;
  					//newborn male individuals
  f,i,
  kids,
  child,
  partner : integer;

begin
with world[x,y] do
  begin //1
  a:=(lambda_0-1)/capacity;		//calculation of crowding-parameter
  if a<0 then a:=0;	//makes sure, that survival does not become
  			//artifically high, if lambda_0 < 1
  sqr_dens:=sqr((fpopsize+mpopsize)/capacity);	//calculation of square-density
  						//for allee-effect
  survival:=(sqr_dens/(sqr_dens+par.allee*par.allee))
  	    /(1+a*(fpopsize+mpopsize));		//survival-probability of
            					//newborns based on theory of
            					//logistic growth


  for f:=1 to fpopsize do    //every female has offspring
    begin //2
    kids:=poisson(lambda*survival);	//number of newborns is poisson-distr.
    					//on the base of logistic
      for child:=1 to kids do	
        begin //3		

		inc(newfpopsize);

        if newmpopsize>MAXM then
          begin
          newmpopsize:=MAXM;
          inc(par.offspr_exceeded);
          end;

        if newfpopsize>MAXF then
          begin
          newfpopsize:=MAXF;
          inc(par.offspr_exceeded);
          end;

        female[f].genetics(newfemale[newfpopsize]);

        end; //3

    end;  //2

    fpopsize:=newfpopsize;		//old females die, young will live
    for i:=1 to fpopsize do
      begin
      female[i]:=newfemale[i];		//old females get overwritten by young
      end;

  end; //1
end;

//______________________________________________________________________________
//-------------------------------------------------------------Reproduction Loop

procedure reproduction_loop;
var
  x,y : integer;
  lambda : real;			//per capita growth rate, calculated by
  					//loggauss distribution
  newfpopsize,
  newmpopsize : integer;		//numbers of newborn individuals

begin

  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
    begin  //1

    with world[x,y] do
      begin //2
      
      ana.pop_b4_rep[x,y]:=fpopsize+mpopsize;
      
      newfpopsize:=0;
      newmpopsize:=0;
      if lambda_0>0 then lambda:=loggauss(lambda_0,sigma)
      		    else lambda:=0;

      if (fpopsize>0) then reproduction(x,y,
      							   lambda,
      							   newfpopsize,
      							   newmpopsize)
        else
        begin
        fpopsize:=0;
        mpopsize:=0;
        end;

      popsize:=fpopsize+mpopsize;
      ana.pop_after_rep[x,y]:=popsize;

      end; //2


    end;  //1

end;

//______________________________________________________________________________
//------------------------------------------------Calculate Dispersal propensity

function disp_prob(x,y : integer; l1a1,l1a2 : real):real;
begin	

  if par.ddd then 
    begin			//calculate dispersal propensity following
  	if density[x,y]>0 then		//Poethke & Hovestadt (2002)
	  begin

	  if (density[x,y] < ((l1a1+l1a2)/2))
		then
		disp_prob:=0
		else
		disp_prob:=1-((l1a1+l1a2)/2)*(1/density[x,y]);

	  end
	  else disp_prob:=0;
	end else
  begin
  	if ((l1a1+l1a2)/2)>0 then disp_prob:=((l1a1+l1a2)/2)
						 else disp_prob:=0;
  end;
end;

//______________________________________________________________________________
//---------------------------------------------------Find Dispersal Target Patch

procedure find_target_patch(const xsource,ysource: integer);
begin

  repeat xtarget:=random(par.max_x)+1 until xtarget<>xsource;
  repeat ytarget:=random(par.max_y)+1 until ytarget<>ysource;
  
end;

//______________________________________________________________________________
//---------------------------------------------------------------------Dispersal

procedure dispersal(const xs,ys,xt,yt : integer;var ind : integer;sex : string);
var
  popsize_b4_f,popsize_b4_m : integer;
begin
  inc(ana.emigrants[xs,ys]);
  inc(ana.emigrants_all);

  popsize_b4_f:=world[xt,yt].fpopsize;
  popsize_b4_m:=world[xt,yt].mpopsize;

  if sex='male' then
    begin
    
      inc(ana.emigrants_m[xs,ys]);

      if random>((world[xs,ys].disp_mort+world[xt,yt].disp_mort)/2) then
        begin				//the mean of dispersal mortality betw.
        inc(world[xt,yt].mpopsize);	//s and t is calculated for direction
        inc(ana.immigrants[xt,yt]);
        inc(ana.immigrants_m[xt,yt]);
        world[xt,yt].male[world[xt,yt].mpopsize]:=world[xs,ys].male[ind];
        world[xt,yt].male[world[xt,yt].mpopsize].dispersed:=true;
        
	  if ((popsize_b4_m=0) and (popsize_b4_f=0)) then
	    begin
        inc(ana.colonized);
		end;

        end;     //if survival

     world[xs,ys].male[ind]:=world[xs,ys].male[world[xs,ys].mpopsize];
     dec(world[xs,ys].mpopsize);

	 if (world[xs,ys].mpopsize+world[xs,ys].fpopsize)=0 then
		inc(ana.extinct);

        dec(ind);

    end;  //if sex=male

  if sex='female' then
    begin
    
      inc(ana.emigrants_f[xs,ys]);

      if random>((world[xs,ys].disp_mort+world[xt,yt].disp_mort)/2) then
        begin
        inc(world[xt,yt].fpopsize);
        inc(ana.immigrants[xt,yt]);
        inc(ana.immigrants_f[xt,yt]);
        world[xt,yt].female[world[xt,yt].fpopsize]:=world[xs,ys].female[ind];
        world[xt,yt].female[world[xt,yt].fpopsize].dispersed:=true;

	    if ((popsize_b4_f=0) and (popsize_b4_m=0)) then
	    begin
        inc(ana.colonized);
		end;

        end;  //if survival

        world[xs,ys].female[ind]:=world[xs,ys].female[world[xs,ys].fpopsize];
        dec(world[xs,ys].fpopsize);
     
	  if (world[xs,ys].mpopsize+world[xs,ys].fpopsize)=0 then
         inc(ana.extinct);

	dec(ind);

    end; //if sex=female

end;

//______________________________________________________________________________
//----------------------------------------------------------------Dispersal-Loop

procedure dispersal_loop;
var
  f,m : integer;		//running variables for females and males
  d : real;			//dispersal propensity calculated by disp_prob
begin

for xsource:=1 to par.max_x do
for ysource:=1 to par.max_y do
  begin  //1
  if ((world[xsource,ysource].fpopsize+world[xsource,ysource].mpopsize)>0) then
    begin //2

    f:=0;					//dispersal of females
    while f<world[xsource,ysource].fpopsize do
	begin //3
      inc(f);
	  
	if not world[xsource,ysource].female[f].dispersed then
	  begin
		  inc(ana.count_loc1_all);
		  ana.sum_loc1_all:=ana.sum_loc1_all
					+((world[xsource,ysource].female[f].loc1_allel1
					+world[xsource,ysource].female[f].loc1_allel2)/2);
		  inc(ana.sum_ind_b4[xsource,ysource]);
		  inc(ana.sum_ind_b4_all);
	  end;

      d:=disp_prob(xsource,ysource,
				   world[xsource,ysource].female[f].loc1_allel1,
				   world[xsource,ysource].female[f].loc1_allel2);

      if ((random<d) and not (world[xsource,ysource].female[f].dispersed)) then
        begin //4
	find_target_patch(xsource,ysource);

        dispersal(xsource,ysource,xtarget,ytarget,f,'female');
		end; //4

      end; //3

    m:=0;					//dispersal of males
    while m<world[xsource,ysource].mpopsize do
	  begin //5
	  inc(m);

	if not world[xsource,ysource].male[m].dispersed then
	  begin
		  inc(ana.count_loc1_all);
		  ana.sum_loc1_all:=ana.sum_loc1_all
					+((world[xsource,ysource].male[m].loc1_allel1
					+world[xsource,ysource].male[m].loc1_allel2)/2);
		  
		  inc(ana.sum_ind_b4[xsource,ysource]);
		  inc(ana.sum_ind_b4_all);
	  end;

	  d:=disp_prob(xsource,ysource,
								 world[xsource,ysource].male[m].loc1_allel1,
								 world[xsource,ysource].male[m].loc1_allel2);

	  if ((random<d) and not (world[xsource,ysource].male[m].dispersed)) then
		begin //6
	find_target_patch(xsource,ysource);

        dispersal(xsource,ysource,xtarget,ytarget,m,'male');
		end; //6
      end; //5

    end; //2

  end; //1

end;

//______________________________________________________________________________
//--------------------------------------------------------------------Extinction

procedure extinction;
var
  x,y : integer;				//running variables
begin

  for x:=1 to par.max_x do
  for y:=1 to par.max_y do
    begin //1

    if (ana.occupied[x,y] and (random<world[x,y].ext_prob)) then
      begin
      world[x,y].mpopsize:=0;
      world[x,y].fpopsize:=0;
      ana.occupied[x,y]:=False;
      dec(ana.meta_occupied);
      end;

    end; //1

end;

//______________________________________________________________________________
//------------------------------------------------------------------------Census

procedure count_patches_pre;
var
	x,y : integer;
begin
	ana.metapop_pre:=0;
	ana.meta_occupied_pre:=0;
	for x:=1 to par.max_x do
	for y:=1 to par.max_y do
	begin
		world[x,y].occ_pre_b4:=world[x,y].occ_pre;
		ana.metapop_pre:=ana.metapop_pre+world[x,y].fpopsize+world[x,y].mpopsize;
		if ((world[x,y].fpopsize+world[x,y].mpopsize)>0) then 
		begin
			world[x,y].occ_pre:=true;
			inc(ana.meta_occupied_pre);
		end else world[x,y].occ_pre:=false;
	end;
end;

procedure count_patches_post;
var
	x,y : integer;
begin
	ana.metapop_post:=0;
	ana.meta_occupied_post:=0;
	for x:=1 to par.max_x do
	for y:=1 to par.max_y do
	begin
		world[x,y].occ_post_b4:=world[x,y].occ_post;
		ana.metapop_post:=ana.metapop_post+world[x,y].fpopsize+world[x,y].mpopsize;
		if ((world[x,y].fpopsize+world[x,y].mpopsize)>0) then 
		begin
			world[x,y].occ_post:=true;
			inc(ana.meta_occupied_post);
		end else world[x,y].occ_post:=false;
	end;
end;

//______________________________________________________________________________
//-------------------------------------------------------------Running Procedure

procedure simulate;
var
  t,repli,x,y : integer;
  mus,sigs : Array[1..101] of double;
  mu,s : integer;
  shelper,mhelper : double;
begin

  shelper:=+0;		// these variables get the temporary values of sigma
  mhelper:=+0;		// and mu for the simulations (see below)
  for s:=1 to 101 do // 101 values are tested
  begin
  sigs[s]:=shelper;
  mus[s]:=mhelper;
  shelper:=shelper+0.02;	// in these steps
  mhelper:=mhelper+0.01;
  end;
  
  for mu:=1 to 101 do
  begin
  par.disp_mort:=mus[mu];
  
  for s:=1 to 101 do
  begin
  par.sigma:=sigs[s];

    for repli:=1 to par.replications do
      begin

//initialize the population----------------------------------------------

      makepop;
//time-loop--------------------------------------------------------------

      t:=0;
      while t<par.tmax do
        begin
        inc(t);
				
//clear disperal-analysis-variables--------------------------------------

        clear_analysis_disp;
				
//calculate patch densities for density-dependent dispersal--------------

				for x:=1 to par.max_x do    	
				for y:=1 to par.max_y do 			
					begin
					if world[x,y].capacity>0 then
					density[x,y]:=(world[x,y].mpopsize+world[x,y].fpopsize)/
								   world[x,y].capacity else density[x,y]:=0;
					end;
				
//sum up population size b4 the last gen.

for x:=1 to par.max_x do
for y:=1 to par.max_y do
begin
  ana.pop_last_gen[x,y]:=world[x,y].mpopsize+world[x,y].fpopsize;
  ana.pop_last_gen_f[x,y]:=world[x,y].fpopsize;
  ana.pop_last_gen_m[x,y]:=world[x,y].mpopsize;
end;
				
//dispersal, sex and extinctions-----------------------------------------

		count_patches_pre;

        dispersal_loop;

        count_patches_post;
        
        calc_occs;

        reproduction_loop;
        extinction;
				
//analyze at certain moments in time-------------------------------------

        if (t > (par.tmax-500)) then analyze(outfile,t,repli);

        end;

      end;
  
  end;
  end;

end;

begin
  init;
  simulate;
  close(outfile);
  close(parinfile);
end.
