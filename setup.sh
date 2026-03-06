#!/usr/bin/env bash

# Make sure we're up to date with upstream and origin.
git pull --no-edit upstream main
git pull --no-edit


read -r -p "Do you want to convert Quarto to Jupyter? (Only do this when starting the worksheet.) [Y/n] " convert
if [[ "${convert:-y}" =~ ^[Yy]$ ]]; then
  for quarto in *.qmd; do
    WORKSHEET="${quarto%.*}"
    QMD="${WORKSHEET}.qmd"
    IPYNB="${WORKSHEET}.ipynb"

    if [[ -f "${IPYNB}" ]]; then
      echo "${IPYNB} already exists. Make sure you save it before proceeding."
      read -r -p "Do you want to back up your existing notebook? [Y/n] " backup
      if [[ "${backup:-y}" =~ ^[Yy]$ ]]; then
        IPYNB_BACKUP="${WORKSHEET}-$(date '+%Y%m%d-%H%M%S').ipynb"
        echo "Backing up ${IPYNB} to ${IPYNB_BACKUP}"
        cp "${IPYNB}" "${IPYNB_BACKUP}"
      fi
    fi

    quarto render "${QMD}" --to ipynb
    python ../../tools/.clean_ipynb.py "${IPYNB}"
  done
fi
