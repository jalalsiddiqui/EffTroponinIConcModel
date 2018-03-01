//Develop Stoichiometric Matrix
//List of reactions and species
//Reaction 
//1.  TnC + Ca ---> TnC-Ca
//2.  TnC-Ca ---> TnC + Ca
//3.  Tnc-Ca + TnI ---> TnC-Ca-TnI
//4.  TnC-Ca-TnI ---> Tnc-Ca + TnI
//5.  TnC-Ca-TnI ---> TnC-TnI + Ca
//6.  TnC-TnI + Ca ---> TnC-Ca-TnI
//7.  TnC-TnI ---> TnC + TnI
//8.  TnI + actin ---> TnI-Actin
//9.  TnI-actin ---> TnI + actin
//10. EGTA + Ca ----> EGTA-Ca
//11. EGTA-Ca --> EGTA + Ca

//12. TnC + Mg ----> TnC-Mg
//13. TnC-Mg ---> TnC + Mg
//14. EGTA + Mg ---> EGTA-Mg
//15. EGTA-Mg ---> EGTA + Mg

//16. TnC-TnI + Mg ---> TnC-TnI-Mg
//17. TnC-TnI-Mg ---> TnC-TnI + Mg
//18. TnC-TnI-Mg ---> TnC-Mg + TnI


//Species
//1.  TnC
//2. Ca
//3. TnC-Ca
//4. TnI
//5. TnC-Ca-TnI
//6. TnC-TnI
//7. actin
//8. TnI-actin
//9. EGTA
//10. EGTA-Ca

//11. Mg
//12. TnC-Mg
//13. EGTA-Mg

//14. TnC-TnI-Mg

v(5) = 1;
v(14) = 0;

//Stochiometric Matrix: S(Species, Reactions) 



function dataMatrix = extractWTData()
    //extracts wt transient occupancy data into dataMatrix
      [fd,SST,Sheetnames,Sheetpos] = xls_open('TOData.xls');
      [dataMatrix,TextInd] = xls_read(fd, Sheetpos);
      mclose(fd);
      
      
    endfunction
    
function dataMatrix1 = extractMutData()
    //extracts truncated TnI transient occupancy data into dataMatrix
      [fd,SST,Sheetnames,Sheetpos] = xls_open('MutData.xls');
      [dataMatrix,TextInd] = xls_read(fd, Sheetpos);
      dataMatrix1 = [dataMatrix(1:10000,1),dataMatrix(1:10000,2:2:8)]
      
      mclose(fd);
    endfunction
    
    
 function dataMatrix1 = extractTNTData()
     //extracts TnT DK210 transient occupancy data into dataMatrix
      [fd,SST,Sheetnames,Sheetpos] = xls_open('TNTMut.xls');
      [dataMatrix,TextInd] = xls_read(fd, Sheetpos);
      mclose(fd);
      dataMatrix1 = [dataMatrix(1:10000,1),dataMatrix(1:10000,2:2:8)]
    endfunction   



function Overlay(tonumber)
    //plots the transient occupancy data when selected
    if tonumber == 1
        //plots the control data
        data = extractWTData()
        realTime = data(1:10000,1);
        realData = data(1:10000,2:2:8);
        plot(realTime,realData);
    end
    if tonumber == 2
        //plots the truncated TnI
        data = extractMutData()
        realTime = data(1:10000,1);
        realData = data(1:10000,2:1:5);
        plot(realTime,realData);
    end
      if tonumber == 3
          //plots the TnT DK210 data
        data = extractTNTData()
        realTime = data(1:10000,1);
        realData = data(1:10000,2:1:5);
        plot(realTime,realData);
    end 
endfunction

function u = u_reactions(k,x)
//This develops a rate vector for each of the species 
u(1) = k(1)*x(1)*x(2);
u(2) = k(2)*x(3);

u(3) = k(3)*x(3)*x(4);
u(4) = k(4)*x(5);

u(5) = k(5)*x(5);
u(6) = k(6)*x(2)*x(6);

u(7) = k(7)*x(6);

u(8) = k(8)*x(4)*x(7);
u(9) = k(9)*x(8);

u(10) = k(10)*x(2)*x(9);
u(11) = k(11)*x(10);

u(12) = k(12)*x(1)*x(11);
u(13) = k(13)*x(12);

u(14) = k(14)*x(9)*x(11);
u(15) = k(15)*x(13);

u(16) = k(16)*x(6)*x(11);
u(17) = k(17)*x(14);
u(18) = k(18)*x(14);
endfunction

function xstore = SolveMassBalances(S,k,x0,h,T)
 //This function uses the scilab ode solver to solve a solve a series
 //of differential equations given the rate constant vector k, the
 //initial conditions vector x0, the interval h, and the time span T
 //steps is time span over interval
 steps = T/h;
//Set a time vector
timev = 1:steps;
timev = h*timev;
//Calculate a calcium dissociation curve
y =  ode(x0,0,timev,list(du_dx,S,k))
xstore = y';
  
endfunction



function xstore = SolveMassBalancesEquil(S,k,x0,h,T)
   //This function uses the scilab ode solver to solve a solve a series
 //of differential equations given the rate constant vector k, the
 //initial conditions vector x0, the interval h, and the time span T
 //steps is time span over interval
 //Only the final x vector at T is output 
    steps = T/h;
    lenx = length(x0);
    t = (1:steps)*h;
    y =  ode(x0,0,t,list(du_dx,S,k));
    xstore = y(1:lenx,steps);  
endfunction

function sensfinal = pCa(S,k,x0,h,T,egtain,v)
//This function given a rate parameter vector k, initial concentrations 
//vector x0, time span T, and interval h, it calculates the calcium 
//sensitivity curve for these parameters


mprintf("Calculating Calcium Sensitivity\n")
//The initial TnC concentration si stored
//The initial concentration vector is stored
initialtnc = x0(1);
initialx0 = x0;
//Calcium and EGTA levels are set at 0
x0(2) = 0;
x0(9) = 0;


//The following steps allow for equilibriation with Mg
mprintf("Equilibrium with magnesium and other species\n")
magequil = SolveMassBalancesEquil(S,k,x0,h,0.50);

//The rate vector is now the magequil
x0 = magequil;
clear magequil;
//The pCa curves are calculated and printed on the screen.  
//Calcium vectors declared
  Ca = [0.0362 .1 .125 .25 .5 .75 1 1.25 1.5 2.5 5 10 20 40 80 160 320 640 1000]'
  Ca = Ca*10^-6;
//new vectors are declared for free ca and flresults
free_ca = zeros(19);
flresults = zeros(19);

//initial calcium sensitivity curve is calculated
for i = 1:19
//calcium is declared in the initial conditions vector
x0(2) = Ca(i);
//the equilibriated amts are stored in xstore
xstore = SolveMassBalancesEquil(S,k,x0,h,T);

//stores fl and freeca results
flresults(i) = v'*xstore;

free_ca(i) = xstore(2);

end

//resets the initial vector

x0 = initialx0;
x0(2) = 0;
x0(9) = egtain;


//equilibriates with mag
magequil = SolveMassBalancesEquil(S,k,x0,h,0.50);

//resets x0 with mg equilibriation
x0 = magequil;


clear magequil;

 
  //declaring the free_ca1 and flresults1 vector
  //these are indexed with 1
free_ca1 = zeros(19);
flresults1 = zeros(19);

//doing the ca sensivitity calculations
for i = 1:19

x0(2) = Ca(i);
xstore = SolveMassBalancesEquil(S,k,x0,h,T);
    
//Selects for fl. if model is two-state vs three-four-state
//k3 = 0 means model is two state
//else it is three or four state    
//stores fl results
flresults1(i) = v'*xstore;
//final free calcium for each run is stored
free_ca1(i) = xstore(2);

end


 pCaStore = -log10(free_ca);

flresultsfinal = flresults1/initialtnc*100



//The final data points pCaStore and flresults1 are output
//A descriptor is given
mprintf("The Calcium Sensitivity plot with an effective tni of %f uM and an egta amount of %f uM with kinetic parameters\n",x0(4)/1e-6, egtain/1e-6)


    mprintf("k1:  %f\n",k(1));
    mprintf("k2:  %f\n",k(2));
    mprintf("k3:  %f\n",k(3));
    mprintf("k4:  %f\n",k(4));
    mprintf("k5:  %f\n",k(5));
    mprintf("k6:  %f\n",k(6));
    mprintf("k7:  %f\n",k(7));
    mprintf("k12:  %f\n",k(12));
    mprintf("k13:  %f\n",k(13));
    mprintf("tni:  %f\n",initialx0(4));
    mprintf("mg:  %f\n",initialx0(11));
    //The curve is output
    mprintf("pCa,  Saturation\n");
    for i = 1:19
        mprintf("%f,  %f\n",pCaStore(i),flresultsfinal(i));
    end
    //The calcium sensitivity pCa50 is calculated
    pCa50 = interp1(flresultsfinal,pCaStore,0.5*max(flresultsfinal));
    mprintf("pCa:  %f\n",pCa50);
    
    //the sensitivity and saturation points are stored in sensfinal
    sensfinal(1) = pCa50;
    sensfinal(2) = max(flresultsfinal)


endfunction



function rodecay = decayRate(S,k,x0,h,T,v,filename)
//given initial conditions x0, rate constants k, time span T and interval h
//the function generates a decay curve and calculates a decay rate based on 
//t1/2
//stores dissociation curve and results in decay_rate.txt
mprintf("Calculating Rate of Dissociation of Calcium with 10 mM EGTA\n")
//store effective TnI concentration
tnieff = x0(4);
//set egta and calcium concentration to 0
x0(9) = 0;
x0(2) = 0;

mprintf("Equilibrium with magnesium and other species\n")
//m0 is the equilibriated concentrations when mg is added
m0 = SolveMassBalancesEquil(S,k,x0,h,0.50);
//200 uM of calcium is added
m0(2) = 200e-6;
//newConc is equilibriated concentrations when 200 uM Ca2+ is added
newConc = SolveMassBalancesEquil(S,k,m0,h,0.50);
//add 10 mM EGTA to the resulting newConc mixture
newConc(9) = 10e-3;
//calculates decay vector xstore
//selects fl results
xstore = SolveMassBalances(S,k,newConc,h,T)*v;    


//declaration of the time vector for the fl. results
steps = T/h;
timev = 1:steps;
timev = h*timev; 

//Calculations of exponential decay
//This is the exponential rate calculator
//the maximum of the curve is determined from the newConc results
//depends on fl. vector
    maxdecay = v'*newConc;


//max decay is displayed
disp('maxdecay');
disp(maxdecay/1e-8)
//the minimum fl. value from the xstore vector is used
mindecay = 0;
//the difference between the the maximum and minimum fl results is calculated
difference_decay = maxdecay - mindecay;
//the time of the half-way point is determined- t1/2
point_half = interp1(xstore,timev,mindecay + 0.5*difference_decay);
//a is the decay rate calculated from t1/2 or point_half
a = log(2)/(point_half);
//write results to decay_rate.txt file
    fid = mopen(filename,"w");
     mfprintf(fid,"Given kinetic parameters:  \n")
    mfprintf(fid,"k1:  %f\n",k(1));
    mfprintf(fid,"k2:  %f\n",k(2));
    mfprintf(fid,"k3:  %f\n",k(3));
    mfprintf(fid,"k4:  %f\n",k(4));
    mfprintf(fid,"k5:  %f\n",k(5));
    mfprintf(fid,"k6:  %f\n",k(6));
    mfprintf(fid,"k7:  %f\n",k(7));
    mfprintf(fid,"TNI:  %f uM\n",tnieff/1e-6);
    
    mfprintf(fid,"The rate of dissociation of calcium is %f s-1\n", a);
    
    mfprintf(fid,"Time(s)            TnC-TnI-Ca\n")
    for i = 1:steps
        
        mfprintf(fid,"%f,  %f\n",i*h,100*xstore(i)/1e-6);
        
    end
    mclose(fid)
    //end file transfer
    
    //output results to screen
    mprintf("Given kinetic parameters:  \n")


    mprintf("k1:  %f\n",k(1));
    mprintf("k2:  %f\n",k(2));
    mprintf("k3:  %f\n",k(3));
    mprintf("k4:  %f\n",k(4));
    mprintf("k5:  %f\n",k(5));
    mprintf("k6:  %f\n",k(6));
    mprintf("k7:  %f\n",k(7));
    mprintf("TNI:  %f uM\n",tnieff/1e-6);
    
    mprintf("The rate of dissociation of calcium is %f s-1\n", a);
    
    //store rate of decay in rodecay
    rodecay = a;

//disp(xstore/1e-6);

endfunction


function to = RunTransOccupancy(S,k,x0,h,T,ca,v)
  //calculates a transient occupancy curve TnC and TnC complexes given a set of
  //kinetic parameters k, initial conditions x0, time span T, interval h, and
  //calcium inputs ca
  //outputs a transient occupancy curve to 'to'
  mprintf("Running Transient Occupancy with 600 uM EGTA and %f uM calcium\n",ca/1e-6)
  //Initial calcium concentration at 0
  x0(2) = 0;
  //EGTA concentration at 600 uM
  x0(9) = 600e-6;
  //Equilibriate with magnesium and store results onto x0
  mprintf("Equilibriating with Magnesium and other species\n")
  x0 = SolveMassBalancesEquil(S,k,x0,h,0.50);
  
  //Calculate a transient curve with 2000 uM calcium to serve as 100% point
  //use maximum as 100% point
  //set calcium to 2000 uM
  x0(2) = 2000e-6;
  //calculate the curve
  //selects fl results
  t2000 = SolveMassBalances(S,k,x0,h,T)*v;  
  //the max at 2000 uM is stored
  max2000 = max(t2000);
  
  //calculate the calcium transient curve based on the calcium input:  ca
  //set x0(2) or calcium levels of x0 to input ca <ca>
  x0(2) = ca;
  //calculate the curve and outputs fl. vector
      to = SolveMassBalances(S,k,x0,h,T)*v/max2000*100;
 
  //determines and displays the amplitude in the curve
  maxto = max(to);
  //displays the transient occupancy amplitude
  mprintf("Transient occupancy amplitude is %f \n",maxto)

endfunction

function kapp = runD(S,k,x0,h,T,tni,egtain,v,tonumber)
    //Runs a series of simulations given rate constants k, initial conditions x0,
    //effective tni concentrations, egta for pCa, and tonumber
    //stores data summary in kapp:  1. saturation 2. pCa50 3. decay rate based
    //on t1/2
    
    //sets effective TnI concentration. This step allows for ease of calculating
    //tni curves
    x0(4) = tni;
    //calculates the pCa curve and displays pCa points on screen
    pc = pCa(S,k,x0,h,2, egtain,v);
    //calculates teh decay curve and displays t1/2 rate on screen
    //outputs results to decay_rate.txt
    dr = decayRate(S,k,x0,h,T,v,'random.txt');
    //stores data summary in kapp:  1. saturation 2. pCa50 3. decay rate based
    //on t1/2
    kapp(1) = pc(2);
    kapp(2) = pc(1);
    kapp(3) = dr;
    
    to12 = RunTransOccupancy(S,k,x0,h/2,T,12.5e-6);
    to25 = RunTransOccupancy(S,k,x0,h/2,T,25e-6);
    to50 = RunTransOccupancy(S,k,x0,h/2,T,50e-6);
    to1000 = RunTransOccupancy(S,k,x0,h/2,T,1000e-6);
    
    fid = mopen('transoccupancy.txt',"w");
    mfprintf(fid,"time, to12, to25, to50, to1000, %f. %f\n", k(1), k(2));
    steps = T/(h/2);
    tprint = 1:steps;
    tprint = tprint*h/2;
    
    
    figure;
    plot(tprint,to12);
    
    plot(tprint,to25);
    plot(tprint,to50);
    plot(tprint,to1000);
    
    for i = 1:steps
        
      mfprintf(fid, "%f, %f, %f, %f, %f\n", tprint(i),to12(i),to25(i),to50(i),to1000(i));
   end
    mclose(fid);
    


    Overlay(tonumber);
  
endfunction

function kapp = runDTO(S,k,x0,h,T,tni,egtain,tonumber,v)
    //Runs a series of simulations given rate constants k, initial conditions x0,
    //effective tni concentrations, egta for pCa, and tonumber
    //stores data summary in kapp:  1. saturation 2. pCa50 3. decay rate based
    //on t1/2
    
    //sets effective TnI concentration. This step allows for ease of calculating
    //tni curves
    x0(4) = tni;
    //calculates the pCa curve and displays pCa points on screen
    pc = pCa(S,k,x0,h,0.60, egtain,v);
    //calculates teh decay curve and displays t1/2 rate on screen
    //outputs results to decay_rate.txt
    dr = decayRate(S,k,x0,h,T,v,random);
    //stores data summary in kapp:  1. saturation 2. pCa50 3. decay rate based
    //on t1/2
    kapp(1) = pc(2);
    kapp(2) = pc(1);
    kapp(3) = dr;
    //Runs a series of transient occupancy curves and outputs data to 
    //'transoccupancy.txt' and plots data onto figure
    //runs data for ca2+ inputs of 12.5 uM, 25 uM, 50 uM, and 1000 uM
    to12 = RunTransOccupancy(S,k,x0,h/2,T,12.5e-6,v);
    to25 = RunTransOccupancy(S,k,x0,h/2,T,25e-6,v);
    to50 = RunTransOccupancy(S,k,x0,h/2,T,50e-6,v);
    to1000 = RunTransOccupancy(S,k,x0,h/2,T,1000e-6,v);
    //opens file transoccupancy.txt and starts writing to file
    fid = mopen('transoccupancy.txt',"w");
    mfprintf(fid,"time, to12, to25, to50, to1000, %f. %f\n", k(1), k(2));
    //declares a time vector tprint to output onto transoccupancy.txt
    steps = T/(h/2);
    tprint = 1:steps;
    tprint = tprint*h/2;
    
   //outputs vectors onto file using following formatting
   //text file is comma separated
    for i = 1:steps
        
      mfprintf(fid, "%f, %f, %f, %f, %f\n", tprint(i),to12(i),to25(i),to50(i),to1000(i));
   end
   //close file
   mclose(fid);
   //create new figure and plot data points on figure
    figure;
    plot(tprint,to12);
    
    plot(tprint,to25);
    plot(tprint,to50);
    plot(tprint,to1000);
    
   
//overlay a dataset for the transient occupancy curve
//1 means control
//2 means trunc TnI
//3 means TnT DK210    
Overlay(tonumber);
  
endfunction

function riseRate(S,k,x0,h,T,tni,ca,v)
//Calculates the calcium association curve given a set of kinetic constants k,
//effective TnI concentration tni, calcium input ca, initial conditions x0, and
//interval h and span T.  
    
  mprintf("Calculating the Rate of Association\n");
  //set initial calcium to zero
  x0(2) = 0;
  //set egta to 0
  x0(9) = 0;
  //reset tni effective concentration
  x0(4) = tni;
  //equil. with magnesium
  mprintf("Equilibriating with magnesium and other species\n");
  //equil. results with magnesium
  x0 = SolveMassBalancesEquil(S,k,x0,h,0.50);
  //set calcium input to ca
  x0(2) = ca;
  //calculates the rise curve
  rise = SolveMassBalances(S,k,x0,h,T)*v;
  //declare a time-rise set of values
  timerise = h*(1:(T/h));
  //determine the max rise pt
  maxrise = max(rise);
  //estimate the rise rate through the half-way point
  point_half = interp1(rise,timerise,0.5*maxrise);
  //calculate the rate through half-way
  a = log(2)/point_half;
  //output the rise rate results
  mprintf("Given kinetic paraemeters\n")
  mprintf("k1:  %f\n",k(1));
    mprintf("k2:  %f\n",k(2));
    mprintf("k3:  %f\n",k(3));
    mprintf("k4:  %f\n",k(4));
    mprintf("k5:  %f\n",k(5));
    mprintf("k6:  %f\n",k(6));
    mprintf("k7:  %f\n",k(7));
    mprintf("TNI:  %f uM\n",x0(4)/1e-6);
    
    mprintf("Rate of Calcium association:  %f s-1\n",a)
endfunction

function dudt = du_dx(t,x,S,k)
//Calculates the species rate vector
dudt = S*u_reactions(k,x);
endfunction

function output = runCurves(S,k,x0,h,T,tni,k4,v, egta)
k(4) = k4;
k(7) = k4;
x0(4) = tni;
pc = pCaRush(S,k,x0,h*100,30,egta,v);
output(1) = pc(1);
output(2) = decayRateRush(S,k,x0,h,0.50,v);
endfunction

S = zeros(8,9);

//Input values for Stoichiometric Matrix.  For all reactions involving finding k_on using probes and chelators.  Use this matrix.  
S(1,1) = -1;
S(2,1) = -1;
S(3,1) = 1;

S(1,2) = 1;
S(2,2) = 1;
S(3,2) = -1;

S(3,3) = -1;
S(4,3) = 0;
S(5,3) = 1;

S(3,4) = 1;
S(4,4) = 0;
S(5,4) = -1;

S(2,5) = 1;
S(5,5) = -1
S(6,5) = 1;



S(2,6) = -1;
S(5,6) = 1
S(6,6) = -1;

S(1,7) = 1;
S(4,7) = 0;
S(6,7) = -1;

S(4,8) = 0;
S(7,8) = -1;
S(8,8) = 1;

S(4,9) = 0;
S(7,9) = 1;
S(8,9) = -1;

S(9,10) = -1
S(2,10) = -1
S(10,10) = 1;

S(9,11) = 1
S(2,11) = 1
S(10,11) = -1;

//Reaction 12

S(1,12) = -1;
S(11,12) = -1;
S(12,12) = 1 ;

//Reaction 13
S(1,13) = 1;
S(11,13) = 1;
S(12,13) = -1 ;

//Reaction 14

S(9,14) = -1
S(11,14) = -1
S(13,14) = 1

//Reaction 15

S(9,15) =  1
S(11,15) = 1
S(13,15) = -1

//Reaction 16

S(6,16) =  -1
S(11,16) = -1
S(14,16) = 1

//Reaction 17

S(6,17) =  1
S(11,17) = 1
S(14,17) = -1

//Reaction 18

S(14,18) = -1;
S(4,18) = 0;
S(12,18) = 1;


//Input Value for Kinetics Constants
k(1) =  2e8//4e7 //1e7//1e8//4e7 //binding of TnC to Ca
k(2) = 1900 //dissociation of Tnc from Ca
k(3) = .9e8 //binding of TncCa to Tni
k(4) = 110 //dissociation of TncCa and Tni
k(5) = 40//330//55//110 //dissociation of Ca from TnCTni
k(6) = 2e8//4e7; //binding of Ca to TnCTnI
k(7) = 110 //dissociation of TnC and TnI
k(8) = 1e8 //Tni binding to actin
k(9) = 200; //Tni dissociating from actin
k(10) = 1.6317e6; //Ca binding to EGTA
k(11) = 0.685; //EGTA dissocation from Ca
k(12) = 1.8e6; //TnC binding to Mg
k(13) = 3000; //TnC dissociating from Mg
k(14) = 7.6923e4; //EGTA binding to Mg
k(15) = 3000; //EGTA dissociating from Mg

k(16) = 1.8e6;//TnC-TnI binding to Mg
k(17) = 3000;//TnC-TnI dissociating from Mg
k(18) = 110;

//Input values for initial conditions

x0(1) = 1e-6;
x0(2) = 25e-6 //1e-6 //Ca
x0(3) = 0;
x0(4) = 8.0e-6
x0(5) = 0;
x0(6) = 0;
x0(7) = 0;
x0(8) = 0;
x0(9) = 500e-6 //500e-6;
x0(10) = 0;
x0(11) = 3e-3;
x0(12) = 0;
x0(13) = 0;
x0(14) = 0;


h = 1e-5; //interval
T = h*3000;

b1(1) = 0.2e-6;
b1(4) = 6e-6;
b1(2) = 200e-6;
b1(11) = 3e-3;
b1(20) = 0;

b2(15) = 10e-6;
b2(4) = 6e-6;
b2(2) = 200e-6;
b2(11) = 3e-3;
b2(20) = 0;


disp('This program will run simulations for the thin filament model given the desired parameters');

disp('The kinetic constants used are as follows:  ');

disp(k);


for i = 1:14
mprintf('dx(%i)/dt = ',i);
for j = 1:17
if S(i,j) ~= 0
    if S(i,j) > 0
        mprintf('+')
    end
    if S(i,j) < 0
        mprintf('-')
    end
mprintf('k(%i)', j)
for m = 1:14
if S(m,j) == -1
mprintf('*x(%i)',m );

end

end



end
end
mprintf(';\n');
end

function sensfinal = pCaRush(S,k,x0,h,T,egtain,v)
//This function given a rate parameter vector k, initial concentrations 
//vector x0, time span T, and interval h, it calculates the calcium 
//sensitivity curve for these parameters


//The initial TnC concentration si stored
//The initial concentration vector is stored
initialtnc = x0(1);
initialx0 = x0;
//Calcium and EGTA levels are set at 0
x0(2) = 0;
x0(9) = 0;


//The following steps allow for equilibriation with Mg

magequil = SolveMassBalancesEquil(S,k,x0,h,0.50);

//The rate vector is now the magequil
x0 = magequil;
clear magequil;
//The pCa curves are calculated and printed on the screen.  
//Calcium vectors declared
  Ca = [0.0362 .1 .125 .25 .5 .75 1 1.25 1.5 2.5 5 10 20 40 80 160 320 640 1000]'
  Ca = Ca*10^-6;
//new vectors are declared for free ca and flresults
free_ca = zeros(19);
flresults = zeros(19);

//initial calcium sensitivity curve is calculated
for i = 1:19
//calcium is declared in the initial conditions vector
x0(2) = Ca(i);
//the equilibriated amts are stored in xstore
xstore = SolveMassBalancesEquil(S,k,x0,h,T);

//stores fl and freeca results
flresults(i) = v'*xstore;

free_ca(i) = xstore(2);

end

//resets the initial vector

x0 = initialx0;
x0(2) = 0;
x0(9) = egtain;


//equilibriates with mag
magequil = SolveMassBalancesEquil(S,k,x0,h,0.50);

//resets x0 with mg equilibriation
x0 = magequil;


clear magequil;

 
  //declaring the free_ca1 and flresults1 vector
  //these are indexed with 1
free_ca1 = zeros(19);
flresults1 = zeros(19);

//doing the ca sensivitity calculations
for i = 1:19

x0(2) = Ca(i);
xstore = SolveMassBalancesEquil(S,k,x0,h,T);
    
//Selects for fl. if model is two-state vs three-four-state
//k3 = 0 means model is two state
//else it is three or four state    
//stores fl results
flresults1(i) = v'*xstore;
//final free calcium for each run is stored
free_ca1(i) = xstore(2);

end


 pCaStore = -log10(free_ca);

flresultsfinal = flresults1/initialtnc*100



//The final data points pCaStore and flresults1 are output
//A descriptor is given

    //The calcium sensitivity pCa50 is calculated
    pCa50 = interp1(flresultsfinal,pCaStore,0.5*max(flresultsfinal));

    
    //the sensitivity and saturation points are stored in sensfinal
    sensfinal(1) = pCa50;
    sensfinal(2) = max(flresultsfinal)


endfunction



function rodecay = decayRateRush(S,k,x0,h,T,v)
//given initial conditions x0, rate constants k, time span T and interval h
//the function generates a decay curve and calculates a decay rate based on 
//t1/2
//stores dissociation curve and results in decay_rate.txt

//store effective TnI concentration
tnieff = x0(4);
//set egta and calcium concentration to 0
x0(9) = 0;
x0(2) = 0;


//m0 is the equilibriated concentrations when mg is added
m0 = SolveMassBalancesEquil(S,k,x0,h,0.50);
//200 uM of calcium is added
m0(2) = 200e-6;
//newConc is equilibriated concentrations when 200 uM Ca2+ is added
newConc = SolveMassBalancesEquil(S,k,m0,h,0.50);
//add 10 mM EGTA to the resulting newConc mixture
newConc(9) = 10e-3;
//calculates decay vector xstore
//selects fl results
xstore = SolveMassBalances(S,k,newConc,h,T)*v;    

//declaration of the time vector for the fl. results
steps = T/h;
timev = 1:steps;
timev = h*timev; 

for i = 1:steps
    xstore(i) = (xstore(i) - xstore(steps))/(v'*newConc-xstore(steps));
end
counter = 0;
for i = 1:steps
    if xstore(i) > 0.1
        counter = counter + 1;
    end
    
end
newvector = log(xstore(1:counter));
[a,b,sig]=reglin(timev(1:counter),newvector');

rodecay = a;



endfunction
