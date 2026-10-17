ARG BASE_IMAGE
FROM ${BASE_IMAGE}
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends gfortran nodejs && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6
COPY world.py run_experiments.py build_pes.py world-browser.js browser_qa.js /app/
COPY vendor /app/vendor
RUN python build_pes.py && gfortran -O2 -fPIC -shared vendor/scalar_lsth.f90 -o vendor/liblsth.so
ENTRYPOINT ["python","-B","run_experiments.py"]
