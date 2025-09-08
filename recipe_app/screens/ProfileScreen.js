import React, { useState, useEffect } from 'react';
import { View, Text, FlatList, StyleSheet, Button } from 'react-native';
import { auth, db } from '../services/firebaseConfig';
import { collection, query, where, onSnapshot } from 'firebase/firestore';

export default function ProfileScreen({ navigation }) {
  const [liked, setLiked] = useState([]);
  const [created, setCreated] = useState([]);
  const uid = auth.currentUser?.uid;

  useEffect(() => {
    if (!uid) return;
    const likedQ = query(collection(db, 'recipes'), where('likedBy', 'array-contains', uid));
    const createdQ = query(collection(db, 'recipes'), where('authorId', '==', uid));

    const unsubLiked = onSnapshot(likedQ, (snap) => {
      setLiked(snap.docs.map((d) => ({ id: d.id, ...d.data() })));    });
    const unsubCreated = onSnapshot(createdQ, (snap) => {
      setCreated(snap.docs.map((d) => ({ id: d.id, ...d.data() })));  });

    return () => { unsubLiked(); unsubCreated(); };
  }, [uid]);

  return (
    <View style={styles.container}>
      <Text style={styles.header}>Liked Recipes</Text>
      <FlatList
        data={liked}
        keyExtractor={(item) => item.id}
        renderItem={({ item }) => <Text style={styles.item}>{item.title}</Text>}
      />
      <Text style={styles.header}>My Recipes</Text>
      <FlatList
        data={created}
        keyExtractor={(item) => item.id}
        renderItem={({ item }) => <Text style={styles.item}>{item.title}</Text>}
      />
      <Button title="Back" onPress={() => navigation.goBack()} color="#0066cc" />
    </View>
  );
}

const styles = StyleSheet.create({
  container: { flex: 1, padding: 20, backgroundColor: 'white' },
  header: { fontSize: 20, marginVertical: 10, color: '#0066cc' },
  item: { paddingVertical: 4, borderBottomWidth: 1, borderBottomColor: '#eee', color: '#0066cc' }
});
