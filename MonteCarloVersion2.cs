using System; 
using System.Collections.Generic;
using System.Text;
using System.Linq;
namespace MyApp
{
    class Program
    {
        static void Main(string[] args)
        {
            int vanderCorput = 0;
            int a = 2;
            int b = 5;
            int anti = 0;
            Console.WriteLine("Welcome.");
            Console.WriteLine("Please input simulation and option parameters.");
            Console.WriteLine("How many paths would you like to simulate?");
            int M = Convert.ToInt32(Console.ReadLine());

            Console.WriteLine("How many steps would you like to simulate?");
            int N = Convert.ToInt32(Console.ReadLine());

            if(N==1){
                Console.WriteLine("I see you only want to simulate one step. Press 1 to engage Van der Corput.");
                vanderCorput = Convert.ToInt32(Console.ReadLine());
                Console.WriteLine("Enter base a. Please ensure you pick a prime number.");
                a = Convert.ToInt32(Console.ReadLine());
                Console.WriteLine("Enter base b. Please ensure you pick a prime number.");
                b = Convert.ToInt32(Console.ReadLine());
            }

            Console.WriteLine("Is this option a call, or is it a put? Please press '0' for a call, and '1' for a put.");
            int call_or_put = Convert.ToInt32(Console.ReadLine());

            Console.WriteLine("What is the strike of the option?");
            double strike = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("In years, what is the tenure of the option?");
            double T = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("As a decimal, what is the current risk free rate?");
            double r = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("What is the current price of the underlying?");
            double S0 = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("What is the current volatility of the underlying?");
            double sigma = Convert.ToDouble(Console.ReadLine()); 

            if(vanderCorput != 1){
                Console.WriteLine("Engage Antithetic Variance Reduction? Press '1' to activate. Press '0' to skip.");
                anti = Convert.ToInt32(Console.ReadLine()); 
            }
            Console.WriteLine("");
            if(vanderCorput==1){
                double[,] vdcPaths = VanDerCorputPathGenerator(a,b,S0,T,sigma,r);
                Console.WriteLine("The value of the option is:");
                Console.WriteLine(OptionPricer(vdcPaths,strike,r,T,call_or_put));
                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("Its greeks are:");
           
                double[] greeks = VDCGreeksCalculator(a,b,S0,strike,T,sigma,r,call_or_put);
                string[] greeksNames = new string[] {"delta","gamma","vega","theta","rho"};
                for(int i = 0; i<greeksNames.Length;i++)
                {
                    Console.WriteLine( greeksNames[i] + ": " + greeks[i]);
                }
            }else{
                double[][,] sim = PathGenerator(M,N,S0,T,sigma,r,anti);
                Console.WriteLine("The value of the option is:");
                Console.WriteLine(OptionPricer(sim[0],strike,r,T,call_or_put));
                Console.WriteLine("");
                Console.WriteLine("The standard error is:");
                double[] optionPrices = new double[M];
                for(int j=0; j<M ; j++)
                {
                    optionPrices[j] = (OptionPricer(sim[0][j,N-1],strike,call_or_put) + OptionPricer(sim[1][j,N-1],strike,call_or_put))/2;
                }
                Console.WriteLine(StandardErrorCalculator(optionPrices));
                Console.WriteLine("");
                Console.WriteLine("Its greeks are:");
           
                double[] greeks = GreeksCalculator(M,N,S0,strike,T,sigma,r,call_or_put,anti);
                string[] greeksNames = new string[] {"delta","gamma","vega","theta","rho"};
                for(int i = 0; i<greeksNames.Length;i++)
                {
                    Console.WriteLine( greeksNames[i] + ": " + greeks[i]);
                }
            }

        }
        static double[] BoxMuller()
        {
            // makes a draw from the normal distribution via the box muller method
            // RETURNS A PAIR OF NUMBERS

            Random r = new Random();
            double x1 = r.NextDouble();
            double x2 = r.NextDouble();
            double z1 = Math.Sqrt(-2*Math.Log(x1))*Math.Cos(2*Math.PI*x2);
            double z2 = Math.Sqrt(-2*Math.Log(x1))*Math.Sin(2*Math.PI*x2);
            double[] normalDraws = {z1,z2};
            return normalDraws;
        }
        static double[] BoxMuller(double x1, double x2)
        {
            // makes a draw from the normal distribution via the box muller method
            // instead of generating random numbers, user passes in numbers instead.

            double z1 = Math.Sqrt(-2*Math.Log(x1))*Math.Cos(2*Math.PI*x2);
            double z2 = Math.Sqrt(-2*Math.Log(x1))*Math.Sin(2*Math.PI*x2);
            double[] normalDraws = {z1,z2};
            return normalDraws;
        }
        //Simulates stock prices
        static double[][,] PathGenerator(int M, int N, double S0, double T, double sigma, double r, int anti)
        {
        
            double t = new double();
            t = T/N;
            

            // initialize with M rows (M paths) and N columns (N time steps)
            double[,] normalDraws = new double[M, N];

        
            // fill it with draws from N(0,1)
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < M; j++)
                {
                    normalDraws[j,i] = BoxMuller()[1];
                }
            }
            
            //initialize stock prices matrix and set its first column to S0
            double[,] stockPrices = new double[M,N];
            double[,] antiPrices = new double[M,N];
            for(int j=0; j < M ; j++)
            {
                stockPrices[j,0] = S0;
                antiPrices[j,0] = S0;

            }


            //step across columns (through time) accroding to geom brownian motion

            if(anti==1){
                for(int i = 0; i<(N-1); i++){
                    for(int j = 0; j < M; j++){
                        stockPrices[j,i+1] = stockPrices[j,i]*Math.Exp((r-Math.Pow(sigma,2)*.5)*t + sigma*Math.Sqrt(t)*normalDraws[j,i]);
                        antiPrices[j,i+1] = antiPrices[j,i]*Math.Exp((r-Math.Pow(sigma,2)*.5)*t - sigma*Math.Sqrt(t)*normalDraws[j,i]);
                    }
                }
            }else{
                for(int i = 0; i<(N-1); i++){
                    for(int j = 0; j < M; j++){
                        stockPrices[j,i+1] = stockPrices[j,i]*Math.Exp((r-Math.Pow(sigma,2)*.5)*t + sigma*Math.Sqrt(t)*normalDraws[j,i]);
                    }
                }

            }
            double[][,] finalArray = new double[2][,];

            if(anti==1){
                finalArray[0] = stockPrices;
                finalArray[1] = antiPrices;
            }else{
                finalArray[0] = stockPrices;
                finalArray[1] = stockPrices;
            }
            


            return(finalArray);
        }
        
        static double[,] VanDerCorputPathGenerator(int a, int b,double S0, double T, double sigma, double r)
        {   
            
            int[] index = Enumerable.Range(1,10000).ToArray();
            int indexHolder = new int();
            int[] baseAReverse = new int[50];
            int[] baseBReverse = new int[50];            
            double sumA = new double();
            double sumB = new double();
            double[] vdcSequenceA = new double[10000];
            double[] vdcSequenceB = new double[10000];
            double[] normalDraws = new double[20000];
            double[,] stockPrices = new double[20000,2];
            
            // for each number from 1 to 20000
            for(int j = 0; j<10000;j++){
                indexHolder = index[j];
                //to base a number
                for(int i = 0; i<50;i++){
                    baseAReverse[i] = indexHolder % a;
                    indexHolder = indexHolder/a;

                    // running tally of the jth vdc number 
                    sumA += baseAReverse[i]*Math.Pow(Convert.ToDouble(a),-(i+1));
                    
                }
                vdcSequenceA[j] = sumA;
                sumA = 0;
            }

            for(int j = 0; j<10000;j++){
                //to base b number
                indexHolder = index[j];
                for(int i = 0; i<50;i++){
                    baseBReverse[i] = indexHolder % b;
                    indexHolder = indexHolder/b;

                    // running tally of the jth vdc number 
                    sumB += baseBReverse[i]*Math.Pow(Convert.ToDouble(b),-(i+1));
                    
                }
                vdcSequenceB[j] = sumB;
                sumB = 0;
            }
            
            //use box muller and save the values.
            for(int i = 0; i<10000; i++){
                BoxMuller(vdcSequenceA[i],vdcSequenceB[i]).CopyTo(normalDraws,2*i);
            }
           

            //generate paths
            for(int j = 0; j < 20000; j++){
                stockPrices[j,0] = S0;
                stockPrices[j,1] = S0*Math.Exp((r-Math.Pow(sigma,2)*.5)*T + sigma*Math.Sqrt(T)*normalDraws[j]);
            }

            return(stockPrices);                
            
        }


        static double OptionPricer(double[,] stockPrices, double strike, double rate, double T, int call_or_put)
        {
            // initalize terminal prices array
            double[] terminalPrices = new double[stockPrices.GetLength(0)];

            // populate array with last column from stockPrices
            for(int j = 0; j<stockPrices.GetLength(0); j++)
            {
            terminalPrices[j] = stockPrices[j,stockPrices.GetLength(1)-1];
            
            }
            
            // calculate option value for each generated path and discount
            double[] payoff = new double[terminalPrices.Length];

            //assume given option is a call.
            if(call_or_put == 1){
                for(int i = 0; i<terminalPrices.Length;i++ )
                {
                    payoff[i] = Math.Max(strike - terminalPrices[i] ,0);
                    payoff[i] = payoff[i]*Math.Exp(-rate*T);
                }
            }else{
                 for(int i = 0; i<terminalPrices.Length;i++ )
                {
                    payoff[i] = Math.Max(terminalPrices[i]-strike ,0);
                    payoff[i] = payoff[i]*Math.Exp(-rate*T);
                }
            }

            // calculate average and return
            double sum = 0;
            foreach(double element in payoff)
            {
                sum+=element;
            }

            return(sum/(payoff.Length));
        }


       
        static double[] GreeksCalculator(int M, int N, double S0, double strike, double T, double sigma, double r, int call_or_put, int anti)
        {
            double perturbationS0 = S0/5;

            double[][,] regularSim = PathGenerator(M,N,S0,T,sigma,r,anti);
            double[][,] perturbedUpSim = PathGenerator(M,N,S0+perturbationS0,T,sigma,r,anti);
            double[][,] perturbedDownSim = PathGenerator(M,N,S0-perturbationS0,T,sigma,r,anti);

            //calculate delta
            double delta = (OptionPricer(perturbedUpSim[0],strike,r,T,call_or_put) - OptionPricer(perturbedDownSim[0],strike,r,T,call_or_put))/(2*perturbationS0);
            if(anti==1){
                double antiDelta = (OptionPricer(perturbedUpSim[1],strike,r,T,call_or_put) - OptionPricer(perturbedDownSim[1],strike,r,T,call_or_put))/(2*perturbationS0);
                delta = (delta+antiDelta)/2;
            }

            //calculate gamma
            double gamma = (OptionPricer(perturbedUpSim[0],strike,r,T,call_or_put) - 2*OptionPricer(regularSim[0],strike,r,T,call_or_put) + OptionPricer(perturbedDownSim[0],strike,r,T,call_or_put))/(Math.Pow(perturbationS0,2));
            if(anti==1){
            double antiGamma = (OptionPricer(perturbedUpSim[1],strike,r,T,call_or_put) - 2*OptionPricer(regularSim[1],strike,r,T,call_or_put) + OptionPricer(perturbedDownSim[1],strike,r,T,call_or_put))/(Math.Pow(perturbationS0,2));
                gamma = (gamma+antiGamma)/2;
            }

            //calculate vega
            double perturbationSigma = sigma/5;
            double[][,] perturbedSigmaUpSim = PathGenerator(M,N,S0,T,sigma+perturbationSigma,r,anti);
            double[][,] perturbedSimgaDownSim = PathGenerator(M,N,S0,T,sigma-perturbationSigma,r,anti);
            double vega = (OptionPricer(perturbedSigmaUpSim[0],strike,r,T,call_or_put) - OptionPricer(perturbedSimgaDownSim[0],strike,r,T,call_or_put))/(2*perturbationSigma);
            if(anti==1){
                double antiVega =  (OptionPricer(perturbedSigmaUpSim[1],strike,r,T,call_or_put) - OptionPricer(perturbedSimgaDownSim[1],strike,r,T,call_or_put))/(2*perturbationSigma);
                vega = (antiVega+vega)/2;
            }


            //calculate theta   
            double perturbationT = T/10;
            double[][,] perturbedTimeSim = PathGenerator(M,N,S0,T+perturbationT,sigma,r,anti);
            double theta = (OptionPricer(perturbedTimeSim[0],strike,r,T+perturbationT,call_or_put) -  OptionPricer(regularSim[0],strike,r,T, call_or_put))/perturbationT;
            if(anti==1){
                double antiTheta = (OptionPricer(perturbedTimeSim[1],strike,r,T+perturbationT,call_or_put) -  OptionPricer(regularSim[1],strike,r,T, call_or_put))/perturbationT;
                theta = (antiTheta+theta)/2;
            }

            //calculate rho
            double perturbationR = r/5;
            double[][,] perturbedRhoUpSim = PathGenerator(M,N,S0,T,sigma,r+perturbationR,anti);
            double[][,] perturbedRhoDownSim = PathGenerator(M,N,S0,T,sigma,r-perturbationR,anti);
            double rho = (OptionPricer(perturbedRhoUpSim[0],strike,r+perturbationR,T, call_or_put) - OptionPricer(perturbedRhoDownSim[0],strike,r-perturbationR,T, call_or_put))/(2*perturbationR);
            if(anti == 1){
                double antiRho = (OptionPricer(perturbedRhoUpSim[1],strike,r+perturbationR,T, call_or_put) - OptionPricer(perturbedRhoDownSim[1],strike,r-perturbationR,T, call_or_put))/(2*perturbationR);
                rho = (antiRho + rho)/2;
            }

            double[] greeks = new double[]{delta,gamma,vega,theta,rho};
            return(greeks);
        }

        static double[] VDCGreeksCalculator(int a, int b, double S0, double strike, double T, double sigma, double r, int call_or_put)
        {
            double perturbationS0 = S0/5;

            double[,] regularSim = VanDerCorputPathGenerator(a,b,S0,T,sigma,r);
            double[,] perturbedUpSim = VanDerCorputPathGenerator(a,b,S0+perturbationS0,T,sigma,r);
            double[,] perturbedDownSim = VanDerCorputPathGenerator(a,b,S0-perturbationS0,T,sigma,r);

            //calculate delta
            double delta = (OptionPricer(perturbedUpSim,strike,r,T,call_or_put) - OptionPricer(perturbedDownSim,strike,r,T,call_or_put))/(2*perturbationS0);
        

            //calculate gamma
            double gamma = (OptionPricer(perturbedUpSim,strike,r,T,call_or_put) - 2*OptionPricer(regularSim,strike,r,T,call_or_put) + OptionPricer(perturbedDownSim,strike,r,T,call_or_put))/(Math.Pow(perturbationS0,2));
     

            //calculate vega
            double perturbationSigma = sigma/5;
            double[,] perturbedSigmaUpSim = VanDerCorputPathGenerator(a,b,S0,T,sigma+perturbationSigma,r);
            double[,] perturbedSimgaDownSim = VanDerCorputPathGenerator(a,b,S0,T,sigma-perturbationSigma,r);
            double vega = (OptionPricer(perturbedSigmaUpSim,strike,r,T,call_or_put) - OptionPricer(perturbedSimgaDownSim,strike,r,T,call_or_put))/(2*perturbationSigma);

            //calculate theta   
            double perturbationT = T/10;
            double[,] perturbedTimeSim = VanDerCorputPathGenerator(a,b,S0,T+perturbationT,sigma,r);
            double theta = (OptionPricer(perturbedTimeSim,strike,r,T+perturbationT,call_or_put) -  OptionPricer(regularSim,strike,r,T, call_or_put))/perturbationT;


            //calculate rho
            double perturbationR = r/5;
            double[,] perturbedRhoUpSim = VanDerCorputPathGenerator(a,b,S0,T,sigma,r+perturbationR);
            double[,] perturbedRhoDownSim = VanDerCorputPathGenerator(a,b,S0,T,sigma,r-perturbationR);
            double rho = (OptionPricer(perturbedRhoUpSim,strike,r+perturbationR,T, call_or_put) - OptionPricer(perturbedRhoDownSim,strike,r-perturbationR,T, call_or_put))/(2*perturbationR);
      

            double[] greeks = new double[]{delta,gamma,vega,theta,rho};
            return(greeks);
        }


        // similar to option pricer, but only accepts and returns a singular double value.
        static double OptionPricer(double underlying, double strike, int call_or_put)
        {
            double optionPrice =  Math.Max((underlying-strike)*Math.Pow(-1,call_or_put),0);
            return (optionPrice);
        }
        static double StandardErrorCalculator(double[] vec)
        {
            //calculate the average value of the vector
            double sum = 0;
            foreach(double element in vec)
            {
                sum += element;
            }
            double avg = sum/vec.Length;

            //calculate the standard deviation of the vector
            double stde = 0;
            foreach(double element in vec)
            {
                stde += Math.Pow((element - avg),2);
            }
            stde = Math.Sqrt(stde/(vec.Length-1));

            //calculate the standard error of the vector
            stde = stde/(Math.Sqrt(vec.Length));
            
            return(stde);

        }
        //create classes for: options, stocks.
        //options: underlying, price, tenure, strike
        //stock: current price, sigma...
    }


}
