## Testing of Python scripts
Python scripts included here require some packages that can be installed in a virtual environment.  This can be accomplished from the root directory of the repo using the following commands:

```bash
python3 -mvenv venv
. venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Run the console script tests using:

```bash
tox
```
