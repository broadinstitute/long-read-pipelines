FROM jupyter/datascience-notebook

########################################
# Update pip and install some basic packages for plotting and data
# analysis:
RUN pip install --upgrade pip \
    && pip install datapackage \
    && pip install folium \
    && pip install --upgrade seaborn \
    && pip install geopandas \
    && pip install tensorflow \
    && pip install nbresuse

# Enable the plugins for jupyter notebooks:
RUN jupyter serverextension enable --py nbresuse --sys-prefix
RUN jupyter nbextension install --py nbresuse --sys-prefix
RUN jupyter nbextension enable --py nbresuse --sys-prefix

RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager

########################################
# Install python bioinformatics / datascience packages:

# Install UMAP clustering package:
RUN pip install llvmlite --ignore-installed
RUN pip install umap-learn

# Install single-cell analysis packages:
RUN pip install scanpy
RUN pip install leidenalg louvain
RUN pip install umap-learn
RUN pip install anndata 
RUN pip install joblib 

# Pysam to read bam/sam/fastX files:
RUN pip install pysam

# Install OpenCV for image processing algorithms:
RUN pip install opencv-python

########################################
# Install libraries for interactive graph visualization:
#RUN pip install bokeh==1.4.0
RUN pip install bokeh
RUN jupyter labextension install @bokeh/jupyter_bokeh

RUN pip install "holoviews[recommended]"
RUN pip install graphviz

RUN pip install plotly==4.9.0
RUN jupyter labextension install jupyterlab-plotly@4.9.0
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager plotlywidget@4.9.0

# Install PyCytoscape
# (https://blog.jupyter.org/interactive-graph-visualization-in-jupyter-with-ipycytoscape-a8828a54ab63)
RUN pip install ipycytoscape
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-cytoscape

# Install the backend support for the graphviz python project:
USER root
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe
RUN apt-get update
RUN apt-get install -y graphviz
RUN apt-get install -y libgraphviz-dev
RUN apt-get install -y graphviz-dev 

RUN apt-get install -y graphviz libgraphviz-dev graphviz-dev pkg-config

RUN echo "jovyan:jovyan" | chpasswd

USER jovyan
RUN pip install pygraphviz
RUN jupyter labextension install @jupyterlab/toc

########################################
# Install some system utilities that would
# be nice to have:

# Install SSW:
RUN mkdir ssw
RUN wget -O ssw-zip.zip https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/master.zip
RUN unzip ssw-zip.zip
RUN cd Complete-Striped-Smith-Waterman-Library-master/src && make libssw.so && mv libssw.so ../../ssw/. && mv ssw_lib.py ../../ssw/.
RUN touch ssw/__init__.py
RUN mv Complete-Striped-Smith-Waterman-Library-master ssw/.
RUN rm ssw-zip.zip 

# Install parasail:
RUN pip install parasail

# Install abPOA:
USER root
RUN apt-get install -y libz-dev graphviz
USER jovyan
RUN wget https://github.com/yangao07/abPOA/releases/download/v1.1.0/abPOA-v1.1.0.tar.gz \
			&& tar -zxf abPOA-v1.1.0.tar.gz \
			&& cd abPOA-v1.1.0 \
			&& make install_py

########################################
# Install some more python packages so we can export plotly images:
RUN pip install -U kaleido

########################################

# When we run the container, start a server:
CMD ["start-notebook.sh"]

