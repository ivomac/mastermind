"""The entry point for playing the Mastermind game.

It handles user interaction, including input for game parameters and guesses,
and displays feedback based on the player's guesses.
"""

import numpy as np

from .engine import Mastermind


def main():
    """Play the Mastermind game.

    Prompts the user for the length of the sequence and the range of numbers,
    then starts the Mastermind game. The user is asked to guess the sequence
    until they get it right.
    """
    n = int(input("Enter the length of the sequence (n): "))
    k = int(input("Enter the range of numbers (0 to k-1): "))
    game = Mastermind(n, k)

    print(f"Welcome to Mastermind! Try to guess the {n}-digit sequence.")
    print(f"Each digit is between 0 and {k-1}. Enter your guess as space-separated numbers.")
    print("✅ indicates correct numbers in the correct position.")
    print("🔄 indicates correct numbers in the wrong position.")
    print("Type 'q' or 'quit' to exit, 'r' or 'restart' for a new sequence.")

    while True:
        guess = input("> ")
        if guess.lower() in {"q", "quit"}:
            break
        if guess.lower() in {"r", "restart"}:
            game.restart()
            print("A new sequence has been generated. Try again!")
            continue

        try:
            guess_list = np.array((int(num) for num in guess.split()), dtype=np.uint8)
        except ValueError:
            print("Invalid input. Please enter numbers only, or 'q' to quit.")
            continue

        if len(guess_list) != n:
            print(f"Please enter exactly {n} numbers.")
            continue

        correct_position, correct_number = game.evaluate_guess(guess_list)
        print(f"✅: {correct_position}  🔄: {correct_number}")

        if correct_position == n:
            print(f"🎉 You guessed the sequence in {game.tries} tries!")
            play_again = (
                input("Do you want to restart with a new sequence? (y/n): ").strip().lower()
            )
            if play_again == "y":
                game.restart()
                print("A new sequence has been generated. Try again!")
            else:
                break
    return


if __name__ == "__main__":
    main()
