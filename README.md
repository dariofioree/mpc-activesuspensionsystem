# mpc-active-suspensions-system
A Model Predictive Control approach for active suspension control problem of a quarter car model
Model-Predictive-Control

Model Predictive Control of a Quarter Car Suspension System with hard constrains

Objective

To model and estimate the parameters of a suspension system and then use MPC control technique within simulation to tune for best performance against different scenarios.

Method

1 - Collect frequency response of a physical Quarter Car Suspension System by inputting a chirp signal to perturb the model
2 - Model the physical system, including Free-Body Diagrams, Equations of Motion, and State Space
3 - Use system identification to estimate parameters and/or get a starting point for a model of the system with explicit coefficient values
4 - Construct MPC controller
5 - Compare results for different scenarios: Rider comfort, Suspension stroke, tyre deflection
6 - Compare results for different speed scenarios: 30 km/h , 60 km/h
Results

MPC performed far better at 30km/h. 
