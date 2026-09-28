"""Paired report text: both languages use identical calculations and structure."""
from contextvars import ContextVar
from functools import wraps
LANG = ContextVar('report_language', default='en')

def tr(en, zh, **values):
    text = zh if LANG.get() == 'zh' else en
    return text.format(**values) if values else text

def localized(fn):
    @wraps(fn)
    def wrapped(*args, lang=None, **kwargs):
        if lang is None:
            return fn(*args, **kwargs)
        if lang not in ('en', 'zh'):
            raise ValueError('Language must be en or zh')
        token = LANG.set(lang)
        try:
            return fn(*args, **kwargs)
        finally:
            LANG.reset(token)
    return wrapped
