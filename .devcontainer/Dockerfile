# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.155.1/containers/ubuntu/.devcontainer/base.Dockerfile

# [Choice] Ubuntu version: bionic, focal
ARG VARIANT="latest"
FROM  jupyter/scipy-notebook:${VARIANT}

USER $NB_UID

# Install Python 3 packages
RUN conda install --quiet --yes \
        "biopython" \
        "flake8" \
        "black" \
    && conda clean --all -f -y \
    && fix-permissions "${CONDA_DIR}" \
    && fix-permissions "/home/${NB_USER}"
