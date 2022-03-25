# Configuration file for the Sphinx documentation builder.
# using commands from https://gist.github.com/matteoferla/ba72ab12a9e5f690277e2e88551773aa
# modified for readthedocs
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# installed in `.readthedocs.yaml`


# -- Project information -----------------------------------------------------

project = 'framenstein'
copyright = '2022, University of Oxford'
author = 'Matteo Ferla etc.'
github_username = 'matteoferla'
github_repository = 'Fragmenstein'


# -- General configuration ---------------------------------------------------

extensions = [
    'readthedocs_ext.readthedocs',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx_toolbox.more_autodoc',
    'sphinx.ext.autodoc',
    #'sphinx.ext.imgconverter',
]

html_static_path = ['_static']
import os
if not os.path.exists('_static'):
    os.mkdir('_static')
for filename in os.listdir('../../images'):  # causing trouble when inspecting files locally in rst
    os.system(f'cp ../../images/{filename} _static/{filename}')

templates_path = ['_templates']
always_document_param_types = True
typehints_defaults = 'braces'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_theme = 'sphinx_rtd_theme'
todo_include_todos = True

# --- add init ---------------------------------------------------------

def skip(app, what, name, obj, would_skip, options):
    if name in ( '__init__',):
        return False
    return would_skip

def setup(app):
    app.connect('autodoc-skip-member', skip)

# --- add mk files ---------------------------------------------------------

import m2r2  # noqa
import os, re

for filename in os.listdir():
    if filename[:4] == 'doc_' or filename[:5] == 'note_':
        os.remove(filename)

repo_base_path = os.path.abspath("../../")

def convert_write(markdown_filename, srt_filename):
    with open(markdown_filename) as fh:
        markdown_block = fh.read()
    markdown_block = re.sub(r"""href=(['"])[./]*images/""", r'href=\1', markdown_block)

    def fix_md_link(match: re.Match) -> str:
        link = match['link']
        if '../' in link or 'documentation/' in link or 'notes/' in link or 'images/' in link:
            pass
        elif 'documentation' in markdown_filename:  # sibling file
            link = 'doc_' + link
        elif 'notes' in markdown_filename:  # sibling file
            link = 'note_' + link
        link = link.replace('../', '')
        link = re.sub(r'^images/', '_static/', link)
        link = re.sub(r'^documentation/notes/', 'note_', link)
        link = re.sub(r'^documentation/', 'doc_', link)
        link = re.sub(r'^notes/', 'note_', link)
        link = re.sub(r'\.md$', '.html', link)
        return f"[{match['label']}]({link})"

    markdown_block = re.sub(r'\[(?P<label>.*?)\]\((?P<link>.*?)\)', fix_md_link, markdown_block)
    rst_block = m2r2.convert(markdown_block)
    with open(srt_filename, 'w') as fh:
        fh.write(rst_block)

convert_write(os.path.join(repo_base_path, 'README.md'), 'introduction.rst')

new_files = {'documentation': [], 'documentation/notes': []}
definitions = (('documentation', 'doc_', 'discussion'),
               ('documentation/notes', 'note_', 'notes'),
               ('documentation/monster', 'doc_', 'discussion'))
for folder, prefix, pagename in definitions:
    for filename in os.listdir(os.path.join(repo_base_path, folder)):
        path = os.path.join(repo_base_path, folder, filename)
        if os.path.isdir(path) or '.md' not in path or 'sphinx' in path:
            continue
        convert_write(path, prefix+filename.replace('.md', '.rst'))
        new_files[pagename].append(prefix+filename.replace('.md', ''))

for pagename in definitions:
    with open(pagename+'.rst', 'w') as fh:
        fh.write(pagename.capitalize()+'\n' +
                 '='*42 +
                 '\n\n.. toctree::\n   :maxdepth: 4\n   :caption: Contents:\n\n' +
                 '\n'.join(['   '+p for p in new_files[pagename]])
                 )
