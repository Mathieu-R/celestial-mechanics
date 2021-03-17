# celestial-mechanics
Simulating Jupiter and Saturn orbiting around the Sun for 5000 years.

### Running simulations
Clone the project
```bash
$ git clone https://github.com/Mathieu-R/celestial-mechanics
```

Create virtual environment
```bash
$ python3 -m venv <env-name>
$ source env/bin/activate
$ python3 -m pip install --upgrade pip
```

Install required packages
```bash
$ python3 -m pip install numpy matplotlib click astropy tqdm
```

Launch simulations
```bash
$ python3 index.py --type=2body --plot=static3D 
$ python3 index.py --type=2body --plot=animated3D
$ python3 index.py --type=3body --plot=static3D
$ python3 index.py --type=3body --plot=animated3D 
```