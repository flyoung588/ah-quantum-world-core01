ARG BASE_IMAGE
FROM ${BASE_IMAGE}
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/app/core03:/app/core02:/app
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends gfortran nodejs && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 pyscf==2.8.0
COPY world.py world-browser.js build_pes.py /app/
COPY vendor /app/vendor
RUN python build_pes.py && gfortran -O2 -fPIC -shared vendor/scalar_lsth.f90 -o vendor/liblsth.so
COPY core02 /app/core02
COPY core03 /app/core03
ENTRYPOINT ["python","-B","core03/run_core03.py"]
