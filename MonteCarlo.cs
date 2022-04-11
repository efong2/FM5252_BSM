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
            Console.WriteLine("Welcome.");
            Console.WriteLine("Please input simulation and option parameters.");
            Console.WriteLine("How many paths would you like to simulate?");
            int M = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("How many steps would you like to simulate?");
            int N = Convert.ToInt32(Console.ReadLine());

            Console.WriteLine("Is this option a call, or is it a put? Please enter '0' for a call, and '1' for a put.");
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
            
            Console.WriteLine("");

            double[,] sim = PathGenerator(M,N,S0,T,sigma,r);
            Console.WriteLine("The value of the option is:");
            Console.WriteLine(OptionPricer(sim,strike,r,T,call_or_put));
            Console.WriteLine("");
            Console.WriteLine("The standard error is:");
            double[] optionPrices = new double[M];
            for(int j=0; j<M ; j++)
            {
                optionPrices[j] = Math.Max((sim[j,N-1]-strike)*Math.Pow(-1,call_or_put),0);
            }
            Console.WriteLine(StandardErrorCalculator(optionPrices));
            Console.WriteLine("");
            Console.WriteLine("Its greeks are:");
           
            double[] greeks = GreeksCalculator(M,N,S0,strike,T,sigma,r,0);
            string[] greeksNames = new string[] {"delta","gamma","vega","theta","rho"};
            for(int i = 0; i<greeksNames.Length;i++)
            {
                Console.WriteLine( greeksNames[i] + ": " + greeks[i]);
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
       
        //Simulates stock prices
        static double[,] PathGenerator(int M, int N, double S0, double T, double sigma, double r)
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
            
            //initialize stock prices matrix
            double[,] stockPrices = new double[M,N];
            for(int j=0; j < M ; j++)
            {
                stockPrices[j,0] = S0;
            }

            //step across columns (through time) accroding to geom brownian motion
            for(int i = 0; i<(N-1); i++){
                for(int j = 0; j < M; j++){
                    stockPrices[j,i+1] = stockPrices[j,i]*Math.Exp((r-Math.Pow(sigma,2)*.5)*t + sigma*Math.Sqrt(t)*normalDraws[j,i]);
                }
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


        //is it better to use monte carlo simulation to estimate the greeks?
        //or is it better to calculate them directly?
        static double[] GreeksCalculator(int M, int N, double S0, double strike, double T, double sigma, double r, int call_or_put )
        {
            double perturbationS0 = S0/5;

            double[,] regularSim = PathGenerator(M,N,S0,T,sigma,r);
            double[,] perturbedUpSim = PathGenerator(M,N,S0+perturbationS0,T,sigma,r);
            double[,] perturbedDownSim = PathGenerator(M,N,S0-perturbationS0,T,sigma,r);

            //calculate delta
            double delta = (OptionPricer(perturbedUpSim,strike,r,T,call_or_put) - OptionPricer(perturbedDownSim,strike,r,T,call_or_put))/(2*perturbationS0);

            //calculate gamma
            double gamma = (OptionPricer(perturbedUpSim,strike,r,T,call_or_put) - 2*OptionPricer(regularSim,strike,r,T,call_or_put) + OptionPricer(perturbedDownSim,strike,r,T,call_or_put))/(Math.Pow(perturbationS0,2));

            //calculate vega
            double perturbationSigma = sigma/5;
            double[,] perturbedSigmaUpSim = PathGenerator(M,N,S0,T,sigma+perturbationSigma,r);
            double[,] perturbedSimgaDownSim = PathGenerator(M,N,S0,T,sigma-perturbationSigma,r);
            double vega = (OptionPricer(perturbedSigmaUpSim,strike,r,T,call_or_put) - OptionPricer(perturbedSimgaDownSim,strike,r,T,call_or_put))/(2*perturbationSigma);

            //calculate theta   
            double perturbationT = T/5;
            double[,] perturbedTimeSim = PathGenerator(M,N,S0,T+perturbationT,sigma,r);
            double theta = (OptionPricer(perturbedTimeSim,strike,r,T+perturbationT,call_or_put) -  OptionPricer(regularSim,strike,r,T, call_or_put))/perturbationT;

            //calculate rho
            double perturbationR = r/10;
            double[,] perturbedRhoUpSim = PathGenerator(M,N,S0,T,sigma,r+perturbationR);
            double[,] perturbedRhoDownSim = PathGenerator(M,N,S0,T,sigma,r-perturbationR);
            double rho = (OptionPricer(perturbedRhoUpSim,strike,r+perturbationR,T, call_or_put) - OptionPricer(perturbedRhoDownSim,strike,r-perturbationR,T, call_or_put))/(2*perturbationR);
            
            double[] greeks = new double[]{delta,gamma,vega,theta,rho};
            return(greeks);
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
