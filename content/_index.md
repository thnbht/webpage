---
date: "2022-10-24"
sections:
- block: about.avatar
  content:
    text: null
    username: admin
  id: about
- block: portfolio
  content:
    buttons:
    - name: All
      tag: '*'
    - name: Data viz
      tag: Data viz
    - name: Immunology
      tag: Immunology
    default_button_index: 0
    filters:
      folders:
      - project
    title: Projects
  design:
    columns: "1"
    flip_alt_rows: false
    view: showcase
  id: projects
- block: contact
  content:
    form:
      formspree:
        id: null
      netlify:
        captcha: false
      provider: netlify
    title: Contact
  design:
    columns: "2"
  id: contact
title: null
type: landing
---
