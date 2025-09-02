from sklearn.linear_model import LinearRegression

def train_model(df):
    # Assumes df has a 'csi' column
    X = df.index.values.reshape(-1, 1)  # Dummy feature: time/order
    y = df['csi'].values
    model = LinearRegression()
    model.fit(X, y)
    return model
