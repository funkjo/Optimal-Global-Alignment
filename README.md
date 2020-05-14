## Optimal Global Alignment Program

School project for Intro to Bioinformatics  
This program reads genomic sequences from `.fasta` files and returns the optimal global alignment and the optimal global alignment score.


## Installation Instructions

#### 1.) Download this project by clicking the big green button in the top right, then clicking Download ZIP. Then unzip the file.

If you do not have python 3 installed on your machine, click the following link and follow the instructions to download python for your operating system: [Download Python](https://realpython.com/installing-python/)  

#### 2.) After this project has been downloaded and unzipped, navigate to the downloads folder in your terminal or console.  

`cd Downloads`

#### 3.) Then navigate into the project directory  

`cd Optimal-Global-Alignment-master`

#### 4.) Next, create a virtual environment in that directory

***On Mac:***  
`python3 -m venv venv`  

***On Windows:***  
`Virtualenv venv`  

#### 5.) Activate your virtual environment

***On Mac:***  
`source venv/bin/activate`  

***On Windows:***  
`venv\Scripts\activate`  

#### 6.) Next, download all of the project dependencies  

`pip install -r requirements.txt`

#### 7.) Run the application

***On Mac:***  
`python3 app.py`  

***On Windows:***  
`app.py`  


## Usage Instructions

This program will ask the user how many `.fasta` files they will be using.  

`> Number of FASTA files: `  

An input outside of 1 or 2 will result in an error and prompt the user to try again.  

Next the program will ask the user for the path or file name of their FASTA files. If the FASTA files are in the same directory as `app.py`, the user can simply write the name of their files. Otherwise, they should specify a path.

