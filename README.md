# Compile

```{bash}
gcc PROJEKT.c libs/src/nrutil.c -Ilibs/include -I/usr/local/dislin -L/usr/local/dislin -Wl,-rpath,/usr/local/dislin -ldislin -lm -o projekt
```

or 

```{bash}
chmod +x run.sh
```

# Run
```{bash}
./projekt
```

or

```{bash}
./run.sh
```