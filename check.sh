#!/usr/bin/env bash

# Define some useful colors
INFO='\033[0;36m'  # Cyan
PASS='\033[0;32m'  # Green
ERROR='\033[0;31m'  # Red
NC='\033[0m'  # No Color

echo -e "${INFO}---BEGIN CORRECTNESS CHECKS---${NC}"

correctness_test() {
  source_file="$1"
  test_file="${2:-test_${source_file}}"
  test_log="${test_file%.*}.log"
  echo "Checking ${source_file} for correctness..."
  if pytest "${test_file}" > "${test_log}"; then
    echo -e "${PASS}All tests passed! See ${test_log} for details.${NC}"
  else
    echo -e "$(
	cat <<- EOF
	${ERROR}${source_file} did not pass all tests.
	See ${test_log} for more details and fix ${source_file}.${NC}
	EOF
  )"
  fi
}

correctness_test gene_finder.py

echo -e "${INFO}---END CORRECTNESS CHECKS---${NC}\n"

echo -e "${INFO}---BEGIN FORMATTING/STYLE CHECKS---${NC}"

# Check formatting of Python files
check_py_format() {
  source_file="$1"
  log_file="format_${source_file%.*}.log"
  echo "Checking the Git exercise for formatting..."
  if black --diff --check "${source_file}" > "${log_file}" 2>&1; then
    echo -e "${PASS}${source_file} formatted correctly!${NC}"
  else
    echo -e "$(
	cat <<- EOF
	${ERROR}${source_file} needs to be reformatted.
	Run "black ${source_file}" to fix, or see ${log_file} for details
	on how to fix the formatting manually.${NC}
	EOF
  )"
  fi
}

check_py_format gene_finder.py
check_py_format test_gene_finder.py

# Check linting of Python files
check_py_style() {
  source_file="$1"
  log_file="lint_${source_file%.*}.log"
  echo "Checking ${source_file} for code style..."
  if pylint "${source_file}" > "${log_file}"; then
    echo -e "${PASS}${source_file} passed all style checks!${NC}"
  else
    echo -e "$(
	cat <<- EOF
	${ERROR}${source_file} has some style errors.
	See ${log_file} for details on what to fix.${NC}
	EOF
  )"
  fi
}

check_py_style gene_finder.py
check_py_style test_gene_finder.py

# # Check Markdown files for Prettier formatting
# check_md_format() {
#   source_file="$1"
#   echo "Checking ${source_file} for formatting..."
#   if prettier --check "${source_file}" 1>/dev/null 2>/dev/null; then
#     echo -e "${PASS}${source_file} formatted correctly!${NC}"
#   else
#     echo -e "$(
# 	cat <<- EOF
# 	${ERROR}${source_file} needs to be reformatted.
# 	Run "prettier -w ${source_file}" to fix.${NC}
# 	EOF
#   )"
#   fi
# }

check_md_format gene-finder-reflection.md
check_md_format reflection.md

echo -e "${INFO}---END FORMATTING/STYLE CHECKS---${NC}"
echo -e "${INFO}All done!${NC}"

