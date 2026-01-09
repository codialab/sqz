import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    return pd, plt


@app.cell
def _(pd):
    df = pd.read_csv("out19.csv", names=["meta-node", "length", "occurrences", "score"])
    df
    return (df,)


@app.cell
def _(df, plt):
    df.plot.scatter(x="length", y="score", logx=True, s=1, c="#1f77b422")
    plt.show()
    return


@app.cell
def _(pd, plt):
    df2 = pd.read_csv("out19_2.csv", names=["node", "l", "o", "s"])
    df2.plot.scatter(x="l", y="s", logx=True, s=1, c="#1f77b422")
    plt.show()
    return (df2,)


@app.cell
def _(df, plt):
    df.plot.scatter(x="occurrences", y="score", logx=True)
    plt.show()
    return


@app.cell
def _(df2, plt):
    df2.plot.scatter(x="o", y="s", logx=True)
    plt.show()
    return


if __name__ == "__main__":
    app.run()
