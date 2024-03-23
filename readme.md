# PISA 
PISA, for protein annotation using a bunch of AI-based tools to predict protein domains, intrinsic disorder region (IDR), PPI sites, DNA binding sites, RNA binding sites, post-translation modification sites, and metal ions binding sites. Nevertheless, PISA aligned the uploaded structure to four structure databases (i.e. Protein Data Bank, SCOP, ECOD, AlphaFoldDB), to find structure similar proteins. By using protein sequence, 3D viewer PDBe Molstar and our custom JavaScript code, user could explore and visualize the protein of their interested. In addition, a Docker images was built for PISA, user could also easy install PISA on his own local systems.

# Installation
User could use web server on  or install locally with a Docker image.

## Download Docker images
```
Docker pull
```

## Dowload data
```
pull
```

## Running the docker image
```
 docker run -it --name web --add-host dockerhost:172.22.148.150  --gpus all  -v `pwd`/apps:/apps -p 9003:9003 -p 9112:9112 nongbaoting/pisa:latest
``
