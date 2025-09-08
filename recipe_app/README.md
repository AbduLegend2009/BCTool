# Gemini Recipe App (Example)

This folder contains a minimal example of a React Native/Expo app that:

* lets users capture a fridge photo
* sends it to Google's Gemini API to extract ingredients
* queries a Firebase `recipes` collection for matching recipes
* allows users to like, rate, and create recipes
* stores likes and ratings in Firestore
* tracks recipes liked or authored by the user

## Setup

1. Install [Expo](https://docs.expo.dev/get-started/installation/) and [Node.js](https://nodejs.org/).
2. Run `npm install` inside this folder to install dependencies.
3. Provide the following environment variables (e.g., using `.env`):

```
FIREBASE_API_KEY=...
FIREBASE_AUTH_DOMAIN=...
FIREBASE_PROJECT_ID=...
FIREBASE_STORAGE_BUCKET=...
FIREBASE_MESSAGING_SENDER_ID=...
FIREBASE_APP_ID=...
GEMINI_API_KEY=...
```

4. Start the app with `npm start`.

## Notes

* The UI uses a simple blue and white theme.
* `HomeScreen` handles camera capture, ingredient editing, and live recipe querying.
* `ProfileScreen` shows recipes liked or created by the authenticated user.
* `CreateRecipeScreen` lets users add new recipes to Firestore.
* `RecipeDetailScreen` displays full recipe info and allows ratings.

This is a skeleton meant to illustrate the architecture; you can expand styling and error handling as needed.
