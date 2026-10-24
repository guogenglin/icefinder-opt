# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from __future__ import annotations

from pathlib import Path
import subprocess

from .log import ICEfinderError, log


def run_command(command: list[str | Path], capture_output: bool = True) -> subprocess.CompletedProcess[str]:
    '''
    Run a tool with consistent diagnostics and retry native multi-threaded BLASTP crashes once.
    '''
    arguments = [str(argument) for argument in command]
    try:
        return subprocess.run(arguments, capture_output = capture_output, text = True, check = True)
    except subprocess.CalledProcessError as initial_error:
        retry_arguments = single_thread_retry(arguments, initial_error.returncode)
        if retry_arguments is None:
            raise ICEfinderError(
                f'{arguments[0]} failed (exit {initial_error.returncode}).\n'
                f'Command: {arguments!r}\n{initial_error.stderr or ""}'
            ) from initial_error
        diagnostic = (initial_error.stderr or '').strip().splitlines()
        query = arguments[arguments.index('-query') + 1] if '-query' in arguments else 'unknown'
        log(f'BLASTP crashed (exit {initial_error.returncode}); retrying once with one thread. '
            f'Query: {query}. Diagnostic: {diagnostic[-1] if diagnostic else "unavailable"}')
        try:
            # A new process reruns the complete query; partial output is never reused.
            return subprocess.run(retry_arguments, capture_output = capture_output, text = True, check = True)
        except subprocess.CalledProcessError as retry_error:
            raise ICEfinderError(
                f'BLASTP failed after a single-threaded retry (exit {retry_error.returncode}).\n'
                f'Command: {retry_arguments!r}\nInitial diagnostic:\n{initial_error.stderr or ""}\n'
                f'Retry diagnostic:\n{retry_error.stderr or ""}'
            ) from retry_error


def single_thread_retry(command: list[str], returncode: int) -> list[str] | None:
    '''
    Build a one-thread retry only for a multi-threaded BLASTP SIGABRT or SIGSEGV failure.
    '''
    if (Path(command[0]).name not in ('blastp', 'blastp.exe') or returncode not in (-6, -11)
            or '-num_threads' not in command):
        return None
    thread_index = command.index('-num_threads') + 1
    if int(command[thread_index]) <= 1:
        return None
    retry = command.copy()
    retry[thread_index] = '1'
    return retry
