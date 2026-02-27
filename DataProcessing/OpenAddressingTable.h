/*
 * OpenAddressingTable.h
 * Summary:
 * - General linear probing hash table.
 * Important notes:
 * - Find and insert return a pointer to the item.
 * TODO items for future work;
 * Currently uses modulo hashing; can be replaced with a
 * stronger mixing function if needed.
 */

#ifndef OPEN_ADDRESSING_TABLE_H
#define OPEN_ADDRESSING_TABLE_H

#include <string>
#include <vector>

template <typename Key, typename Value>

class OpenAddressingTable {
protected:
    // Enumerated types declaration
    enum state {EMPTY, TAKEN};

    // Struct representing items in the table
    struct Item {
        Key key{};
        Value value{};
        state status = EMPTY;
    };

    // Variables
    unsigned long numItems; // Number of items in the table
    std::vector<Item> items;

    /**
     * @brief Hash key function for mixing.
     *
     * Method to better distribute keys across the table and reduce collisions.
     * Currently, it uses std::hash and modulo, but can be replaced with a stronger mixing function if needed.
     *
     * @param key The key to be hashed.
     * @return The hash key used to determine the index for the key in the table.
     * TODO: Replace with stronger mixing function
     */
    virtual size_t hashKey(const Key& key) const {
        return std::hash<Key>{}(key) % items.size();
    }

    /**
     * @brief Rebuilds and resizes the hash table to accommodate increased requirements.
     *
     * This method is called when the load factor exceeds 0.5, resizing the table
     * and reinserting all prior items in the table.
     */
    virtual void rehash() {
        std::vector<Item> preRehash = items;
        items.clear();
        numItems = 0;
        items.resize(nextPrime(preRehash.size() * 2));

        for (size_t i = 0; i < preRehash.size(); ++i) {
            if (preRehash[i].status == TAKEN) {
                insert(preRehash[i].key).first = preRehash[i].value;
            }
        }
    }

    // HOOK FUNCTIONS FOR KmerTable

    /**
     * @brief Exists to be overridden for kmer increment
     * @param value
     */
    virtual void onDuplicate(Value& value) {
        // Default behavior: do nothing (keep existing value)
    }

    /**
     * @brief Exists to be overridden for initial kmer insertion
     * @param value
     */
    virtual void onInitial(Value& value) {
        // Default behavior: do nothing (keep existing value)
    }

    /**
     * @brief Finds the next prime number greater than or equal to n
     * @param n Starting integer to find the next prime
     * @return Next prime number
     */
    static int nextPrime(int n) {
        auto isPrime = [](int x) {
            if (x < 2) return false;
            for (int i = 2; i*i <= x; ++i)
                if (x % i == 0) return false;
            return true;
        };
        while (!isPrime(n)) ++n;
        return n;
    }

public:

    class iterator {
    private:
        OpenAddressingTable* table;
        size_t index;

        /**
         * @brief Advances through the vector until the next item is found.
         */
        void advanceToNextValid() {
            while (index < table->items.size() &&
                   table->items[index].status != TAKEN) {
                ++index;
            }
        }

    public:
        /**
         * @brief
         * @param inputTable Input table to be iterated through
         * @param index
         */
        iterator(OpenAddressingTable* inputTable, size_t inputIndex)
            : table(inputTable), index(inputIndex) {
            advanceToNextValid();
        }

        /**
         * @brief Overload the ++ operator to advance through the vector.
         *
         * Skips empty items, returns a pointer to the next occupied item.
         *
         * @return
         */
        iterator& operator++() {
            ++index;
            advanceToNextValid();
            return *this;
        }

        /**
         * @brief Used in while loops to terminate once whole table is accessed.
         * @param other The other integer being compared.
         * @return Returns a boolean, true if they are not equal.
         */
        bool operator!=(const iterator& other) const {
            return index != other.index;
        }

        /**
         * @brief Used to access a specific item in the table based on the class index.
         * @return Returns a reference to an item.
         */
        Item& operator*() {
            return table->items[index];
        }
    };

    /**
     * @brief
     * @return
     */
    iterator begin() {
        return iterator(this, 0);
    }

    /**
     * @brief
     * @return Returns an iterator object.
     */
    iterator end() {
        return iterator(this, items.size());
    }

    // OPEN ADDRESSING TABLE METHODS
    virtual ~OpenAddressingTable() = default;

    /**
     * @brief Main constructor for OpenAddressingTable.
     *
     * Initializes the hash table with a specified initial size, which is adjusted to the next prime number to
     * help reduce collisions.
     * The number of items is initialized to zero.
     *
     * @param initialSize initial size for construction of the hash table.
     */
    OpenAddressingTable(int initialSize = 101) {
        items.resize(nextPrime(initialSize));
        numItems = 0;
    }

    // Getter method

    /**
     * @brief Returns the number of items in the table.
     * @return numItems The number of items in the table.
     */
    size_t getNumItems() const {
        return numItems;
    }

    /**
     * @brief Inserts an item into the table with the specified key.
     *
     * Inserts an item into the table using a linear probing strategy to resolve collisions.
     * If the key already exists, the onDuplicate method is called to handle the duplicate key scenario.
     *
     * @param key
     * @return Value& Reference to the value associated with the inserted key.
     */
    virtual std::pair<Value&, bool> insert(const Key& key) {

        // Makes sure references are correct after rehashing
        if (numItems + 1 > items.size()/2) {
            rehash();
        }

        size_t index = hashKey(key);
        size_t i = 0;

        while (items[index].status == TAKEN && items[index].key != key) {
            i++;
            index = (index + i) % items.size(); // linear probing
        }

        if (items[index].status == TAKEN) {
            // Return reference to an existing value
            onDuplicate(items[index].value);
            return {items[index].value, false};
        }
        items[index].key = key;
        items[index].value = Value{}; // default construct
        items[index].status = TAKEN;
        ++numItems;

        onInitial(items[index].value);

        return {items[index].value, true};
    }

    /**
     * @brief Searches the hash table for an entry with the specified key.
     *
     * Performs a lookup using the same probing strategy as insertion (linear probing)
     * starting from the hashed index of the key.
     *
     * @param key Key to search for.
     * @return Pointer to the value associated with key if found; otherwise nullptr.
     */
    virtual const Value* find(const Key& key) const {
        size_t index = hashKey(key);
        size_t i = 0;
        while (items[index].status != EMPTY) {
            if (items[index].status == TAKEN && items[index].key == key)
                return &items[index].value;
            i++;
            index = (index + i) % items.size();
        }
        return nullptr;
    }

};
#endif //OPEN_ADDRESSING_TABLE_H