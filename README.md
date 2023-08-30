# BachelorThesis

MATLAB, Genetic Algorithms, Deep Neural Networks, Profiling, Parallel Computing, IS/OS data, Walk Forward Analysis

Reach out to get the full thesis as PDF. 

## About
This Bachelor thesis evaluates some technical strategies of algorithmic trading and compares those with share price predictions using Deep Neural Networks. The daily data of 20 years of 100 shares of the S&P500 index is used which ensures solid backtesting for each strategy. 

<img width="150" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/a9f0fa89-d287-4758-bae9-9e37d2ab1d5e">

This image describes the methodology. Evaluating the performance of a technical or neural network strategy was done by using input data to optimise parameters or a model. After that, signals were generated, and trades were calculated. For both, the technical and machine learning strategy, IS/OS and the walk forward analysis are used. 


## Result

The result of this thesis is that technical strategies outperformed neural network strategies. However, no strategy could reliably outperform the share the strategy was applied on for all results. The best strategy, Double7, outperformed its underlying benchmark only 41 of 100 times. The average return a year for that strategy was 17% and for the underlying share 22,5%. Thus, for long term investing, the buy and hold strategy is preferable. 

<img width="600" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/26bfc996-7288-4396-b321-167102851cbb">

This image visualises the clear winner in diversified long-term investing - buy and hold. 


## Technical Strategies
The following strategies were used: MACD, Double7, and a simple and complex Bollinger Band Strategy. The parameters were choosen (1.) without optimising and (2.) by optimising using a Genetic Algorithm. 

To visualise, I explain and visualise the Double7 strategy: The condition for a long trade is that the current close price must be greater than the 200-period EMA. A long signal occurs when the current low is lower than the lowest low of the previous 7 periods. After that, the sell long triggers when the current high is higher than the high of the last 7 periods.

<img width="600" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/8278052e-19c3-4ea7-80ae-76c5f7175f97">

Visualisation to explain the Double7 strategy 

<img width="450" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/2acc2212-aaa9-49a9-848f-d7b065b89be0">

Result of Double7 on the M3 Company stock 


## Neural Networks
The fitrnet and fitnet libraries of MATLAB are used for creating Deep Nerual Networks. The activation function, the training function, and number of nodes in the first and second layer are (1.) decided by a rule of thumb; hidden layer nodes should be 2/3 of the input nodes plus the size if the output layer, and (2.) optimised by Genetic Algorithm. RMSE was used as cost function, and Kfold Cross Validation was used to prevent overfitting.  

<img width="600" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/5aa9a5a0-9c1e-459c-b5d0-a8a79e14d470">

The result of optimising the first and the second layer of the DNN is not that clear, but there is a tendency that the first layer prefers less nodes, and the second layer prefers more nodes as seen in the image above. 

<img width="600" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/65b6d74c-be71-40c7-b140-d02194e9b51c">

After making trades with the optimised model, some results can be seen in the image above. 


## Risk Management
For risk management, I used (1.) a fixed Stop Loss (SL) and Take Profit (TP), (2.) include the standard deviation of the previous days of the closing price, and (3.) a trailing SL. The last one is raised when the closing price is higher than before. 

<img width="200" alt="image" src="https://github.com/RobbsX/BachelorThesis/assets/79597633/91af4552-3934-4f30-9fd5-6e25feebe83a">

Trailing SL


## Further Work
However, for some shares, using a strategy could outperform the underlying share consistently. When using Neural Networks, it is important to develop a model that fits to a certain market, that is adapted over time and controlled manually for unexpected risks. Also, using attention-based LSTMs might be a potential improvement to predict stock data solely with its price. By adding other independent variables, I expect the result to improve too. 


## Disclaimer
This is not a investing advice but only shows the content of the thesis. 


## Special Thanks
This project was done with ~20h per week of work over 9 months and one hour of review session with my mentor. I improved my skills in coding, maths, machine learning and trading. 
Thank you for this immersive time, Gismondi R.! 
