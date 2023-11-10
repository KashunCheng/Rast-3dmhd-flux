import pickle


def save_user_defined_function_effects(d: dict, filename: str):
    with open(f'{filename}.ue', 'wb') as f:
        pickle.dump(d, f)


def load_user_defined_function_effects(old: dict, filename: str) -> dict:
    with open(f'{filename}.ue', 'rb') as f:
        return {**pickle.load(f), **old}
