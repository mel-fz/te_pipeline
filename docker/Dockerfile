FROM mambaorg/micromamba:1.5.8

USER root
ENV MAMBA_ROOT_PREFIX=/opt/conda
COPY stellarscope_env.yml star_env.yml /tmp/

RUN micromamba create -y -n stellarscope -f /tmp/stellarscope_env.yml \
 && micromamba create -y -n star -f /tmp/star_env.yml \
 && micromamba clean -a -y \
 && rm /tmp/*.yml

ENV PATH=/opt/conda/envs/stellarscope/bin:/opt/conda/envs/star/bin:$PATH
ENTRYPOINT []

# fail the build if any tool is missing
RUN STAR-avx2 --version && stellarscope --version && samtools --version | head -1 \
 && python -c "import pandas, scipy, pysam; print('python libs ok')"
