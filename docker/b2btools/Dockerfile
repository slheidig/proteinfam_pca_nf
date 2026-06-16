FROM python:3.12-slim

LABEL org.opencontainers.image.title="b2btools"
LABEL org.opencontainers.image.version="3.0.8"
LABEL org.opencontainers.image.description="b2bTools + pandas + matplotlib + scikit-learn for og-b2bpca pipeline"
LABEL org.opencontainers.image.source="https://pypi.org/project/b2bTools/"

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    b2bTools==3.0.8 \
    pandas==2.2.3 \
    matplotlib==3.9.4 \
    scikit-learn==1.6.1 \
    numpy==1.26.4

RUN useradd --create-home --shell /bin/bash b2buser
USER b2buser
WORKDIR /home/b2buser

CMD ["b2bTools", "--help"]
