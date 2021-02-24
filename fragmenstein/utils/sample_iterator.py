'''
https://stackoverflow.com/questions/12581437/python-random-sample-with-a-generator-iterable-iterator
'''
import random


def iter_sample_fast(iterable, samplesize):
    results = []
    iterator = iter(iterable)
    # Fill in the first samplesize elements:
    try:
        for _ in range(samplesize):
            results.append(next(iterator))
    except StopIteration:
        raise ValueError("Sample larger than population.")
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterator, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items
    return results
